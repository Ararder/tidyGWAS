

randomizer <- function(end=10) {
  sample(c(1:10),nrow(tbl), replace=TRUE) == 1
}
mock_dbsnp <- function(){
  # set .env to parent.frame() so the mockings persist in the parent function
  load(test_path("fixtures/b38.rds"))
  load(test_path("fixtures/b37.rds"))
  local_mocked_bindings(
    map_to_dbsnp = function(..., build) {
      if(build == 37) return(b37)
      if(build == 38) return(b38)
    }, .env = parent.frame())

  local_mocked_bindings(
    infer_build = function(...) {
      37
    },
    .env = parent.frame()
    )

  local_mocked_bindings(
    get_bsgenome = function(...) {
      list()
    },
    .env = parent.frame()
  )
}


