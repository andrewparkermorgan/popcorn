# s3_dispatch.R
# Hook up S3 dispatch for generics that may not be defined elsewhere

glance <- function(...) UseMethod("glance")