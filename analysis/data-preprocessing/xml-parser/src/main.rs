use quick_xml::events::Event;
use quick_xml::Reader;
use regex::Regex;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::BufReader;

// Function to identify genes from the definition string
fn identify_genes_from_definition(definition: &str) -> String {
    let definition_lower = definition.to_lowercase();
    let mut identified_genes = Vec::new();

    // S segment genes
    if definition_lower.contains("glycoprotein") || definition_lower.contains("gpc") {
        identified_genes.push("GPC");
    }
    if definition_lower.contains("nucleoprotein") || definition_lower.contains(" np ") ||
       definition_lower.contains("(np)") || definition_lower.contains("np gene") ||
       definition_lower.contains("nucleocapsid") {
        identified_genes.push("NP");
    }

    // L segment genes
    if definition_lower.contains("zinc") || definition_lower.contains("ring") ||
       definition_lower.contains("zinc") || definition_lower.contains("(z)") ||
       definition_lower.contains(" z ") || definition_lower.contains("z gene") {
        identified_genes.push("Z");
    }
    if definition_lower.contains("polymerase") || definition_lower.contains(" l ") ||
       definition_lower.contains("l gene") || definition_lower.contains("rna polymerase") {
        identified_genes.push("L");
    }

    if identified_genes.is_empty() {
        "NA".to_string()
    } else {
        identified_genes.join("-")
    }
}

#[derive(Debug, Default, Clone)]
struct GenBankRecord {
    accession: String,
    sequence: String,
    definition: String,
    genes_from_definition: String,
    latitude: Option<f64>,
    longitude: Option<f64>,
    lat_lon_raw: Option<String>,
    organism: Option<String>,
    strain: Option<String>,
    host: Option<String>,
    geo_loc_name: Option<String>,
    collection_date: Option<String>,
    segment: Option<String>,
    country: Option<String>,
    isolate: Option<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.xml> <output.csv>", args[0]);
        std::process::exit(1);
    }

    let input_file = &args[1];
    let output_file = &args[2];

    eprintln!("Reading XML file: {}", input_file);
    let file = File::open(input_file)?;
    let buf_reader = BufReader::new(file);
    let mut reader = Reader::from_reader(buf_reader);
    reader.config_mut().trim_text(true);

    let mut records = Vec::new();
    let mut current_record = GenBankRecord::default();
    let mut current_element = String::new();
    let mut current_qualifier_name = String::new();
    let mut qualifiers: HashMap<String, String> = HashMap::new();
    let mut in_sequence = false;

    // Regex for parsing coordinates with cardinal directions
    let coord_regex = Regex::new(r"(\d+\.?\d*)\s*([NSEW])")?;

    let mut buf = Vec::new();
    loop {
        match reader.read_event_into(&mut buf) {
            Err(e) => {
                eprintln!("Error at position {}: {:?}", reader.buffer_position(), e);
                break;
            }
            Ok(Event::Eof) => break,

            Ok(Event::Start(ref e)) => {
                let name = e.name();
                let name = std::str::from_utf8(name.as_ref())?;
                current_element = name.to_string();

                if name == "INSDSeq" {
                    in_sequence = true;
                    current_record = GenBankRecord::default();
                    qualifiers.clear();
                }
            }

            Ok(Event::End(ref e)) => {
                let name = e.name();
                let name = std::str::from_utf8(name.as_ref())?;

                if name == "INSDSeq" {
                    in_sequence = false;

                    // Process lat_lon coordinates
                    let lat_lon_candidates = ["lat_lon", "latitude_longitude", "lat_long", "coordinates"];
                    for candidate in &lat_lon_candidates {
                        if let Some(lat_lon_str) = qualifiers.get(*candidate) {
                            current_record.lat_lon_raw = Some(lat_lon_str.clone());

                            // Extract coordinates with cardinal directions
                            let matches: Vec<_> = coord_regex.captures_iter(lat_lon_str).collect();

                            if matches.len() >= 2 {
                                // Parse latitude (first coordinate)
                                if let (Some(lat_val), Some(lat_dir)) = (matches[0].get(1), matches[0].get(2)) {
                                    if let Ok(mut lat) = lat_val.as_str().parse::<f64>() {
                                        // South latitudes are negative
                                        if lat_dir.as_str() == "S" {
                                            lat = -lat;
                                        }
                                        current_record.latitude = Some(lat);
                                    }
                                }

                                // Parse longitude (second coordinate)
                                if let (Some(lon_val), Some(lon_dir)) = (matches[1].get(1), matches[1].get(2)) {
                                    if let Ok(mut lon) = lon_val.as_str().parse::<f64>() {
                                        // West longitudes are negative
                                        if lon_dir.as_str() == "W" {
                                            lon = -lon;
                                        }
                                        current_record.longitude = Some(lon);
                                    }
                                }
                            }
                            break;
                        }
                    }

                    // Extract other qualifiers
                    current_record.organism = qualifiers.get("organism").cloned();
                    current_record.strain = qualifiers.get("strain").cloned();
                    current_record.host = qualifiers.get("host").cloned();
                    current_record.geo_loc_name = qualifiers.get("geo_loc_name").cloned();
                    current_record.collection_date = qualifiers.get("collection_date").cloned();
                    current_record.segment = qualifiers.get("segment").cloned();
                    current_record.country = qualifiers.get("country").cloned();
                    current_record.isolate = qualifiers.get("isolate").cloned();

                    records.push(current_record.clone());
                }
                current_element.clear();
            }

            Ok(Event::Text(e)) => {
                if !in_sequence {
                    continue;
                }

                let text = std::str::from_utf8(&e)?;
                let text = text.trim().to_string();
                if text.is_empty() {
                    continue;
                }

                match current_element.as_str() {
                    "INSDSeq_primary-accession" => {
                        current_record.accession = text;
                    }
                    "INSDSeq_definition" => {
                        current_record.definition = text.clone();
                        current_record.genes_from_definition = identify_genes_from_definition(&text);
                    }
                    "INSDSeq_sequence" => {
                        current_record.sequence = text;
                    }
                    "INSDQualifier_name" => {
                        current_qualifier_name = text;
                    }
                    "INSDQualifier_value" => {
                        if !current_qualifier_name.is_empty() {
                            qualifiers.insert(current_qualifier_name.clone(), text);
                            current_qualifier_name.clear();
                        }
                    }
                    _ => {}
                }
            }

            _ => {}
        }
        buf.clear();
    }

    eprintln!("Parsed {} sequences", records.len());

    // Write to CSV
    eprintln!("Writing to CSV: {}", output_file);
    let mut writer = csv::Writer::from_path(output_file)?;

    // Write header
    writer.write_record(&[
        "accession",
        "sequence",
        "definition",
        "genes_from_definition",
        "latitude",
        "longitude",
        "lat_lon_raw",
        "organism",
        "strain",
        "host",
        "geo_loc_name",
        "collection_date",
        "segment",
        "country",
        "isolate",
    ])?;

    // Write data
    for record in records {
        writer.write_record(&[
            &record.accession,
            &record.sequence,
            &record.definition,
            &record.genes_from_definition,
            &record.latitude.map_or(String::new(), |v| v.to_string()),
            &record.longitude.map_or(String::new(), |v| v.to_string()),
            &record.lat_lon_raw.unwrap_or_default(),
            &record.organism.unwrap_or_default(),
            &record.strain.unwrap_or_default(),
            &record.host.unwrap_or_default(),
            &record.geo_loc_name.unwrap_or_default(),
            &record.collection_date.unwrap_or_default(),
            &record.segment.unwrap_or_default(),
            &record.country.unwrap_or_default(),
            &record.isolate.unwrap_or_default(),
        ])?;
    }

    writer.flush()?;
    eprintln!("Done!");

    Ok(())
}