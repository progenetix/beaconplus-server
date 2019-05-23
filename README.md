## BeaconResponse - a Perl CGI for Beacon+ queries on GA4GH compatible databases

### Query parameters

* datasetId
   * the dataset to be queried
   * example: `arraymap`
   * default: `arraymap`
* assemblyId
   * the genome assembly
   * example: `GRCh38`

### Database naming

The script uses some naming conventions for databases and collections:

* MongoDB database (in our implementation)
   * datasetId
* collections
   * `individuals`
   * `biosamples`
   * `callsets`
   * `variants`

### Example use, web call:

SNV query on dataset "dipg"

* http://beacon.arraymap.org/beacon/beaconplus-server/beaconresponse.cgi?datasetId=dipg&referenceName=17&assemblyId=GRCh38&start=7733516&referenceBases=G&alternateBases=A&filters=icdot-C71.7

CNV query (defaults do dataset "arraymap") with bio-metadata component ("filters")

* http://beacon.progenetix.org/beacon/beaconplus-server/beaconresponse.cgi?datasetId=arraymap&referenceName=9&assemblyId=GRCh38&startMin=19,500,000&startMax=21,975,098&endMin=21,967,753&endMax=24,500,000&alternateBases=DEL&filters=ncit:C3224
