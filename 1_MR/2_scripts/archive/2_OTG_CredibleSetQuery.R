# Title: Open Target Genetics GWASCredibleSets Query
# Description: This script defines a function to query the Open Target Genetics
#              GraphQL API for GWAS credible sets associated with a given variant.
# Date: 2025-09-18

# --- Installation ---
# If you haven't installed the required packages, uncomment and run the following lines:
# install.packages("ghql")
# install.packages("jsonlite")
# install.packages("httr")

# --- Load Libraries ---
library(ghql)
library(jsonlite)
library(httr)

#' Convert rsID to Open Targets variantID format (chr_pos_ref_alt)
#'
#' This helper function uses the Ensembl REST API to look up an rsID and
#' convert it to the chromosome_position_reference_alternate format.
#'
#' @param rsId A string representing the rsID (e.g., "rs2476601").
#' @return A string in the format "chr_pos_ref_alt" or NULL if conversion fails.
convert_rsid_to_variant_id <- function(rsId) {
  ensembl_url <- paste0("https://rest.ensembl.org/variation/human/", rsId)
  cat("Contacting Ensembl API to resolve rsID:", rsId, "\n")

  response <- tryCatch(
    {
      GET(url = ensembl_url, add_headers(Accept = "application/json"))
    },
    error = function(e) {
      message("Error contacting Ensembl API: ", e$message)
      return(NULL)
    }
  )

  if (is.null(response) || http_status(response)$category != "Success") {
    message(
      "Failed to retrieve data from Ensembl for ",
      rsId,
      ". Status: ",
      http_status(response)$reason
    )
    return(NULL)
  }

  data <- content(response, "parsed")

  # Find the primary mapping (usually the first one on a canonical chromosome)
  primary_mapping <- NULL
  for (mapping in data$mappings) {
    if (grepl("^\\d+$|^X$|^Y$", mapping$seq_region_name)) {
      primary_mapping <- mapping
      break
    }
  }

  if (is.null(primary_mapping)) {
    message("Could not find a valid chromosomal mapping for ", rsId)
    return(NULL)
  }

  # Allele string is usually in "REF/ALT" format
  alleles <- strsplit(primary_mapping$allele_string, "/")[[1]]
  if (length(alleles) < 2) {
    message("Could not parse reference and alternate alleles for ", rsId)
    return(NULL)
  }

  # Construct the variant ID
  variant_id <- paste(
    primary_mapping$seq_region_name,
    primary_mapping$start,
    alleles[1],
    alleles[2],
    sep = "_"
  )

  return(variant_id)
}


#' Fetch GWAS Credible Sets from Open Target Genetics
#'
#' This function queries the Open Target Genetics GraphQL API to retrieve
#' GWAS credible set data for a specific genetic variant.
#'
#' @param variantId A string representing the variant ID (e.g., "1_10177_A_AC").
#' @param size An integer specifying the number of results to return per page. Defaults to 10.
#' @param index An integer for the page index (0-based). Defaults to 0.
#'
#' @return A list containing the parsed JSON response from the API.
#'         Returns NULL if the query fails.
#' @examples
#' \dontrun{
#'   credible_sets <- get_gwas_credible_sets(variantId = "1_10177_A_AC", size = 5)
#'   # To view the structure of the result
#'   str(credible_sets)
#'   # To view as pretty JSON
#'   cat(jsonlite::toJSON(credible_sets, auto_unbox = TRUE, pretty = TRUE))
#' }
get_gwas_credible_sets <- function(variantId, size = 1000, index = 0) {
  # --- New: Handle rsID conversion ---
  if (startsWith(variantId, "rs")) {
    cat(
      "Input is an rsID. Attempting to convert to chr_pos_ref_alt format...\n"
    )
    converted_id <- convert_rsid_to_variant_id(variantId)
    if (is.null(converted_id)) {
      message("Failed to convert rsID '", variantId, "'. Aborting query.")
      return(NULL)
    }
    cat("Successfully converted '", variantId, "' to '", converted_id, "'\n")
    variantId <- converted_id # Overwrite with the converted ID
  }

  # --- 1. Setup GraphQL Client ---
  # Define the API endpoint
  api_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  # Create a new GraphQL client
  cli <- GraphqlClient$new(url = api_url)

  # --- 2. Define the GraphQL Query ---
  # Create a new query object
  qry <- Query$new()

  # Define the query string. This is the exact query from your request.
  # The query is named 'GWASCredibleSetsQuery' for clarity.
  query_string <- '
    query GWASCredibleSetsQuery($variantId: String!, $size: Int!, $index: Int!) {
      variant(variantId: $variantId) {
        id
        referenceAllele
        alternateAllele
        gwasCredibleSets: credibleSets(studyTypes: [gwas], page: { size: $size, index: $index }) {
          count
          rows {
            studyLocusId
            pValueMantissa
            pValueExponent
            beta
            finemappingMethod
            confidence
            variant {
              id
              chromosome
              position
              referenceAllele
              alternateAllele
            }
            study {
              traitFromSource
              id
              diseases {
                name
                id
                therapeuticAreas {
                  name
                  id
                }
              }
            }
            locus(variantIds: [$variantId]) {
              rows {
                posteriorProbability
              }
            }
            locusSize: locus {
              count
            }
            l2GPredictions(page: { size: 1, index: 0 }) {
              rows {
                target {
                  id
                  approvedSymbol
                }
                score
              }
            }
          }
        }
      }
    }
  '
  qry$query('GWASCredibleSetsQuery', query_string)

  # --- 3. Define Query Variables ---
  # These variables will be passed to the GraphQL query
  variables <- list(
    variantId = variantId,
    size = size,
    index = index
  )

  # --- 4. Execute the Query ---
  cat("Querying Open Targets API for variant:", variantId, "\n")

  # The result is returned as a raw JSON string, so we need to parse it
  response <- tryCatch(
    {
      cli$exec(qry$queries$GWASCredibleSetsQuery, variables)
    },
    error = function(e) {
      message("An error occurred during the API request: ", e$message)
      return(NULL)
    }
  )

  if (is.null(response)) {
    return(NULL)
  }

  # --- 5. Process and Return the Result ---
  # Parse the JSON string into an R list
  result_data <- jsonlite::fromJSON(response, flatten = TRUE)

  cat("Query successful.\n")
  return(result_data$data)
}

# # --- New Example with rsID ---
# cat("\n\n--- Testing rsID Conversion ---\n")
# # This rsID (rs2476601) is for a well-known SNP in the PTPN22 gene, associated with Type 1 Diabetes and Rheumatoid Arthritis.
# example_rs_id <- "rs2476601"
# credible_set_data_from_rsid <- get_gwas_credible_sets(
#   variantId = example_rs_id,
#   size = 5,
#   index = 0
# )

# # Check if the query returned data and print a summary
# if (!is.null(credible_set_data_from_rsid)) {
#   # Print the entire structure of the returned object
#   cat("\n--- Result Structure from rsID ---\n")
#   str(credible_set_data_from_rsid, max.level = 3)

#   if (
#     !is.null(credible_set_data_from_rsid$variant$gwasCredibleSets$rows) &&
#       nrow(credible_set_data_from_rsid$variant$gwasCredibleSets$rows) > 0
#   ) {
#     cat("\n--- Sample Extracted Data from rsID ---\n")
#     first_row_rs <- credible_set_data_from_rsid$variant$gwasCredibleSets$rows[
#       1,
#     ]
#     trait_rs <- first_row_rs$study.traitFromSource

#     cat(paste("Example Trait for", example_rs_id, ":", trait_rs, "\n"))
#   }
# }
