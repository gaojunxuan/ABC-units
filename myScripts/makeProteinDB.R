myDB <- dbAddProtein(myDB,
                     jsonlite::fromJSON("./myScripts/MBP1_MIXOS.json"))
myDB <- dbAddTaxonomy(myDB,
                      jsonlite::fromJSON("./myScripts/MIXOStaxonomy.json"))
