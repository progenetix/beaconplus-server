## BeaconResponse - a Perl CGI for Beacon+ queries on GA4GH compatible databases

### Query parameters

* datasetId
   * the dataset to be queried
   * example: `arraymap`, which will expand to `arraymap_ga4gh` as database name
   * default: `arraymap`
* assemblyId
   * the genome assembly
   * example: `GRCh38`

### Database naming

The script uses some naming conventions for databases and collections:

* MongoDB database (in our implementation)
   * datasetId`_ga4gh`
* collections
   * `individuals`
   * `biosamples`
   * `callsets`
   * `variants`

### Example use, web call:

SNV query on dataset "dipg" with phenotype response:

* http://progenetix.org/beacon/query/?eferenceName=17&referenceBases=G&variants.alternate_bases=A&variants.start=7577121&dataset_ids=dipg&phenotypes=1

CNV query on dataset "dipg", with phenotype response:

* http://progenetix.org/beacon/query/?eferenceName=9&alternateBases=DEL&startMin=19000000&startMax=21984490&endMin=21984490&endMax=25000000&dataset_ids=dipg&phenotypes=1

CNV query (defaults do dataset "arraymap") with bio-metadata component and phenotype response:

* http://progenetix.org/beacon/query/?eferenceName=chr9&alternateBases=DEL&startMin=20000000&startMax=21984490&endMin=21984490&endMax=25000000&biosamples.bio_characteristics.ontology_terms.term_id=NCIT:C3058&biosamples.bio_characteristics.ontology_terms.term_id=NCIT:C3059&phenotypes=1
