args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
d <- read.dcf(file.path(root,"DESCRIPTION"))[1,]
problems <- character()
if(grepl("maintainer@example.org|First.*Last|your-org",paste(d,collapse=" ")))
  problems <- c(problems,"Replace placeholder contact metadata with a usable project mailbox.")
if(!identical(unname(d["Config/HeteroSpaSim/MetadataReady"]),"true"))
  problems <- c(problems,"MetadataReady must be true only after human verification.")
if(length(problems)) {
  writeLines(c("PUBLIC_RELEASE_NOT_READY",problems))
  quit(status=2L)
}
writeLines("PUBLIC_METADATA_GATE_PASS (mailbox ownership is human-confirmed, not network-verified)")
