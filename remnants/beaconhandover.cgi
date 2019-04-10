#!/usr/bin/perl

# Beacon+ support scripts
# © 2017 Michael Baudis: m@baud.is

use strict;
use CGI::Carp qw(fatalsToBrowser);
use CGI qw(:standard param *table);

use JSON;
use MongoDB;
use MongoDB::MongoClient;
use Data::Dumper;
use UUID::Tiny;

=pod

Script for linking callset ids from Beacon+ query to the original data

Please see the associated beaconresponse.md

=cut

our $tempdb     =   'progenetix';
our $tmpcoll    =   'querybuffer';

#if (! -t STDIN) { print 'Content-type: application/json'."\n\n" }

# parameters
our $access_id  =   param('accessid');
our $cgi        =   new CGI;


_add_nice_form();

#_add_form();


################################################################################
# subs #########################################################################
################################################################################

sub _add_nice_form {

=pod

Copied & modified html page from the beacon-ui repository.

=cut

  print <<END;
Content-type: text/html


<!DOCTYPE html>
<html lang="en">

<head>

  <meta charset="UTF-8">
  <title>Beacon+ - Feature evaluation for the GA4GH Beacon project</title>
  <meta name="google-site-verification" content="D_DC_Zxy3cCmrst5uaScIDD7daZ2_HO7Vhq2ZjfSIMo" />
  <meta name="google-site-verification" content="is7-_LmgwVxd7M4_XIn6psLPjlZyk58QuP-DN7FC44g" />
  <link rel="stylesheet" href="/beaconplus-ui/css/bootstrap/css/bootstrap.min.css">
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
  <link rel="stylesheet" href="/beaconplus-ui/css/style.css">
  <script src="https://code.jquery.com/jquery-2.1.4.min.js"></script>
  <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>

</head>

<body>


<div class="container">

  <div class="jumbotron" style="margin-bottom: 10px; padding:  20px 20px 30px 20px;">
    <img src="/beaconplus-ui/img/beacon_logo_25x70.png" style="float: right;" />
    <h2 style="margin: 0px;">Beacon<span style="vertical-align: super; color: red; font-weight: 800;">+</span></h2>
  </div>

  <div class="panel panel-default" >
    <div class="panel-body">
      <div class="alert alert-info" style="min-height: 40px; padding-top: 10px;">
        <p>
          This is an implementation of a Beacon "handover" concept, in which a
          Beacon query response additionally delivers an "accessid" value. This
          value represents an pointer to an internal representation of the query
          resulte (i.e. callsets, biosamples, metadata ...), which can then
          be accessed after authentication.
          The "handover" scenario separates the standard qualitative ("yes"|"no")
          or quantitative ("n matches") Beacon response from a data delivery
          mechanism.
        </p>
        <p>
          The current implementation exemplifies some possible scenarios:
          <ul>
            <li>
              providing a histogram of regional gain/loss frequencies (DUP, DEL)
              for samples with structural variation data
            </li>
            <li>
              returning data of the associated callsets which matched the Beacon
              query (this is for feature demonstration only...)
            </li>
            <li>
              returning the metadata (diagnoses etc.) of the biosamples
              from which the matching callsets were derived
            </li>
          </ul>
        </p>
        <p>
          This demonstrator does not implement authentication procedures yet;
          login & password fields can be left empty.
        </p>
      </div>
    </div>
  </div>

  <div class="panel panel-default" style="min-height: 260px">

    <div class="panel-body">

      <form class="form-horizontal" id="handover-form" method="post" action="/beaconplus-server/beacondeliver.cgi">

        <input type="hidden" name="accessid" id="accessid" value="$access_id" />

        <div class="form-group">
            <label for="todo" class="col-sm-2 control-label">Handover Action</span></label>
            <div class="col-sm-4">
              <select id="todo" name="do" class="form-control">
                <option value="histogram" selected="selected">Plot DUP/DEL histogram</option>
                <option value="callsets">Export Callset Data</option>
                <option value="biosamples">Export Biosample Data</option>
                <option value="variants">Export Variants Data</option>
              </select>
            </div>
            <input type="checkbox" name="jsonpretty" value="1" id="jsonpretty"> pretty JSON
        </div>

        <div class="form-group">
          <label for="login" class="col-sm-2 control-label">Login</label>
          <div class="col-sm-4">
            <input type="text" class="form-control" name="login" id="login" placeholder="login name" />
          </div>
        </div>

        <div class="form-group">
          <label for="pass" class="col-sm-2 control-label">Password</label>
          <div class="col-sm-4">
            <input type="password" class="form-control" name="pass" id="pass" placeholder="password" />
          </div>
        </div>

        <div class="form-group">
            <div class="col-sm-offset-2">
                <div class="col-sm-2">
                  <button class="btn btn-info" target="_blank">Process Data</button>
                </div>
            </div>
        </div>

      </form>

    </div>
  </div>

  <div class="panel panel-default" style="border: none;">
    <div class="panel-body">
      <div style="display: flex; justify-content: space-between;">
        <a href="http://arraymap.org" target="_blank"><img src="/beaconplus-ui/img/arraymap_black_180.png" style="height: 25px; margin: 5px;" /></a>
        <a href="http://progenetix.org" target="_blank"><img src="/beaconplus-ui/img/progenetix_black_180.png" style="height: 25px; margin: 0px 5px 5px 5px;" /></a>
        <div style="min-width: 300px;">
          This Beacon implementation is developed by the <a href="http://wiki.progenetix.org/Wiki/BaudisgroupIMLS/" target="_blank">Computational Oncogenomics Group</a> at the <a href="http://www.imls.uzh.ch/baudis/" target="_blank">University of Zurich</a>, with support from the <a href="https://www.sib.swiss/sibt/" target="_blank">SIB Technology group</a> and <a href="https://www.elixir-europe.org" target="_blank">ELIXIR</a>.
        </div>
        <a href="http://www.uzh.ch" target="_blank"><img src="/beaconplus-ui/img/uzh_logo_pos_160x50.png" style="height: 25px; margin: 5px;" /></a>
        <a href="http://elixir-europe.org" target="_blank"><img src="/beaconplus-ui/img/elixir_logo_grey_65x50.png" style="height: 25px; margin: 5px;" /></a>
        <a href="http://sib.swiss" target="_blank"><img src="/beaconplus-ui/img/sib_logo_65x50.png" style="height: 25px; margin: 5px;" /></a>
      </div>
    </div>
  </div>

</div>

<!-- StatCounter -->
<script type="text/javascript">
  var sc_project=3556072;
  var sc_invisible=1;
  var sc_security="018d1bf6";
  var scJsHost = (("https:" == document.location.protocol) ?
  "https://secure." : "http://www.");
  document.write("<sc"+"ript type='text/javascript' src='" +
  scJsHost + "statcounter.com/counter/counter.js'></"+"script>");
</script>
<noscript>
  <div class="statcounter">
    <a title="Web Analytics"
    href="http://statcounter.com/" target="_blank"><img
    class="statcounter"
    src="//c.statcounter.com/3556072/0/018d1bf6/1/" alt="Web
    Analytics"></a>
</div>
</noscript>
<!-- End of StatCounter -->

<!-- Google -->
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
  ga('create', 'UA-572981-2', 'auto');
  ga('send', 'pageview');
</script>
<!-- End of Google -->


</body>
</html>

END



}

################################################################################

1;
