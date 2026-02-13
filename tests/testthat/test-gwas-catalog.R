test_that("multiplication works", {
  skip()

  false_id <- "GCST99999999"
  from_gwas_catalog(false_id)
  from_gwas_catalog("GCST90441306")
  from_gwas_catalog("GCST90468178")
  from_gwas_catalog("GCST90245992")
  from_gwas_catalog("GCST90451106")


})

test_that(".scrape keeps directories and filters by pattern", {
  html <- paste0(
    "<html><body>",
    "<a href=\"harmonised/\">harmonised</a>",
    "<a href=\"results.tsv\">results</a>",
    "<a href=\"notes.pdf\">notes</a>",
    "<a href=\"../parent\">parent</a>",
    "<a href=\"?query\">query</a>",
    "</body></html>"
  )

  mock_fetch <- function(url) list(content = charToRaw(html))

  result <- with_mocked_bindings(
    tidyGWAS:::`.scrape`("http://example.com", pattern = "\\.(tsv|yaml)$"),
    curl::curl_fetch_memory = mock_fetch
  )

  expect_equal(result, c("harmonised/", "results.tsv"))
})

test_that("scrape_dir retains harmonised files and applies patterns", {
  base_links <- c("harmonised/", "main.tsv", "main.tbi", "readme.txt")
  harmonised_links <- c("harmonised.tsv.gz", "meta.yaml", "other.log")

  mock_scrape <- function(url, pattern) {
    if (grepl("/harmonised/?$", url)) {
      return(harmonised_links)
    }
    base_links
  }

  files <- with_mocked_bindings(
    tidyGWAS:::scrape_dir("http://example.com"),
    tidyGWAS:::`.scrape` = mock_scrape
  )

  expect_equal(
    files,
    c(
      "http://example.com/main.tsv",
      "http://example.com/main.tbi",
      "http://example.com/readme.txt",
      "http://example.com/harmonised/harmonised.tsv.gz",
      "http://example.com/harmonised/meta.yaml",
      "http://example.com/harmonised/other.log"
    )
  )
})
