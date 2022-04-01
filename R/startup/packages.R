library(nat)
library(natverse)
library(neuprintr)
library(nat.templatebrains)
library(nat.jrcbrains)
library(fafbseg)
library(elmr)
library(neuronbridger)
library(lubridate)
library(filesstrings)

register_saalfeldlab_registrations()

#establish neuprint database
conn = neuprint_login(server= "https://neuprint.janelia.org/",
                      token= "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6ImtlbGxvMjJlQG10aG9seW9rZS5lZHUiLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGgzLmdvb2dsZXVzZXJjb250ZW50LmNvbS9hL0FBVFhBSnoxcWZuYXZJTHRHd1FyWGY2OEZGYVdTcFdhU09IMTIwTlRsSFU9czk2LWM_c3o9NTA_c3o9NTAiLCJleHAiOjE4MTA1MTA0MDd9.7IWz4boFJ2iDUmr0zNIjFXGNhSa1KWESateQovBExEE")
options(neuprint_dataset="hemibrain:v1.2.1")