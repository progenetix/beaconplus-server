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

* http://beacon.arraymap.org/beacon/beaconplus-server/beaconresponse.cgi?datasetId=dipg&referenceName=17&assemblyId=GRCh38&start=7733516&referenceBases=G&alternateBases=A&biosamples.biocharacteristics.type.id=pgx:icdot:c71.7&phenotypes=1

CNV query (defaults do dataset "arraymap") with bio-metadata component and phenotype response:

* http://beacon.progenetix.org/beacon/beaconplus-server/beaconresponse.cgi?datasetId=arraymap&referenceName=9&assemblyId=GRCh38&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&alternateBases=DEL&biosamples.biosamples.biocharacteristics.type.id=ncit:c3224&phenotypes=1
