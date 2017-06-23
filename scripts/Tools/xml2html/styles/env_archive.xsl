<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" indent="yes" 
    doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"
    doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN" />

  <xsl:template match="//file/components">
    <html>
      <head>
	<title>CESM CASEROOT env_archive.xml file</title>
	<style type="text/css">
	  body {
          font-family:'Open Sans', Arial, sans-serif;
          font-size:14px;
          font-weight:300;
          line-height:1.6em;
          color:#656565;
          width: 90%;
          align: center;
	  }

	  .container {
          padding: 0 30px 0 30px;
          padding-top: 0px;
          padding-right: 30px;
          padding-bottom: 30px;
          padding-left: 30px;
          position: relative;
	  }
      
	  dt { color:#656565; font-weight: bold }

	  .accordionItem h2 { margin: 0; font-size: 1.1em; padding: 0.4em; color: #000; background-color: #E3E4E6; border-bottom: 1px solid #66d; width: 90%}
	  .accordionItem h2:hover { cursor: pointer; }
	  .accordionItem div { margin: 0; padding: 1em 0.4em; background-color: #eef; border-bottom: 1px solid #66d; width: 90%}
	  .accordionItem.hide h2 { color: #000; background-color: #E3E4E6; width: 90%}
      
	  .detailsBtn { background-color: #E3E4E6; }

	  #variable { color:#66d; font-size: 1.2em; font-weight: bold; display: inline-block; margin-right: 30px; width: 550px; }
	  #emphasis { color:#656565; font-weight: bold; display: inline-block; }
	  #small { color:#9f9f9f; display: inline-block; }

	</style>

	<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js"></script>

	<script type="text/javascript">
	  <xsl:text disable-output-escaping="yes">
	//<![CDATA[ 
	var accordionItems = new Array();

	function init() {
	  // Grab the accordion items from the page
	  var divs = document.getElementsByTagName( 'div' );
	  for ( var i = 0; i < divs.length; i++ ) {
	    if ( divs[i].className == 'accordionItem' ) accordionItems.push( divs[i] );
	  }
	  // Assign onclick events to the accordion item headings
	  for ( var i = 0; i < accordionItems.length; i++ ) {
	    var h2 = getFirstChildWithTagName( accordionItems[i], 'H2' );
	    h2.onclick = toggleItem;
          }
          // Hide all accordion item bodies 
          for ( var i = 0; i < accordionItems.length; i++ ) {
            accordionItems[i].className = 'accordionItem hide';
            $(accordionItems[i]).find('div').slideUp();
          }
        }

        function toggleItem() {
          var itemClass = this.parentNode.className;
          // Hide all items
          for ( var i = 0; i < accordionItems.length; i++ ) {
            accordionItems[i].className = 'accordionItem hide';
            $(accordionItems[i]).find('div').slideUp();
          }
          // Show this item if it was previously hidden
          if ( itemClass == 'accordionItem hide' ) {
            this.parentNode.className = 'accordionItem';
            $(this).parent().find('div').slideDown();
          }
        }

        function getFirstChildWithTagName( element, tagName ) {
          for ( var i = 0; i < element.childNodes.length; i++ ) {
            if ( element.childNodes[i].nodeName == tagName ) return element.childNodes[i];
          }
        }

        function hideAll() {
        // Hide all accordian items
          for ( var i = 0; i < accordionItems.length; i++ ) {
            accordionItems[i].className = 'accordionItem hide';
            $(accordionItems[i]).find('div').slideUp();
          }
        }

        function showAll() {
        // show all accordian items
          for ( var i = 0; i < accordionItems.length; i++ ) {
            accordionItems[i].className = 'accordionItem hide';
            $(accordionItems[i]).find('div').slideDown();
          }
        }
        //]]>
	</xsl:text>
	</script>
      </head>

      <body onload="init()">
	<div class="container">
	  <h1>CESM CASEROOT env_archive.xml file</h1>
	  <p style="font-size: 0.9em;">
	    Model Version: CESM2.0<br/>
	    HTML created on: 2017-06-21
	  </p>

	  <h3>Description</h3>
	  <p>
	    The short-term archiver, case.st_archive, reads these settings in order to<br/>
	    determine the rules for migrating files out of the RUNDIR and into the DOUT_S_ROOT<br/>
	    location while preserving a complete set of necessary restart files in the RUNDIR.
	  </p>
	  <p>
	    These rules are specified in the CASEROOT/env_archive.xml file and rely on python<br/>
	    regular expression filename matches to determine where files should be moved or copied.<br/>
	    Please see <a href="http://www.cesm.ucar.edu/models/cesm2.0/cesm/filename_conventions_cesm.html">
	      CESM2 Output Filename Conventions</a> for details regarding filenames. 
	  </p>
	  <p>
	    The short-term archiver also opens a component restart NetCDF file looking for the<br/>
	    value of restart history variable name &lt;rest_history_varname&gt; to determine<br/>
	    the last component history file required for restarts that must remain as a copy<br/>
	    in the RUNDIR. 
	  </p>

	  <h3>Interface Tools</h3>
	  The env_archive.xml can only be modified manually. Users should run the system <i>xmllint</i><br/>
	  to ensure the validity of the XML schema as follows:
	  <pre>
	    xmllint --schema $CIMEROOT/config/xml_schemas/config_archive.xsd env_archive.xml
	  </pre>

	  <h3>Navigation Instuctions</h3>
	  <p>
	  Clicking on the component name or class will display additional descriptive information.<br/>  
	  Click on the "Show Details" button and then cntl+F key to search for specific strings in this file.
	  </p>

	  <button class="detailsBtn" onclick="showAll();">Show Details</button>
	  <button class="detailsBtn" onclick="hideAll();">Hide Details</button>

	  <br/><br/>
	  <div class="accordionItem">
	    <h2><span id="variable">File Description Header &lt;header&gt;</span></h2>
	    <div>
	      <pre>
		<xsl:value-of select="//file/header" />
	      </pre>
	    </div>
	  </div>
	    
	  <xsl:apply-templates />
	</div>
      </body>
    </html>
  </xsl:template>

  <xsl:template match="comp_archive_spec">
    <div class="accordionItem">
      <h2>Component Archive Specification &lt;comp_archive_spec&gt; :<br/>
	Component Class &lt;compclass&gt; : <span id="variable"><b><xsl:value-of select="@compclass" /></b></span><br/>
	Component Name &lt;compname&gt; : <span id="variable"><b><xsl:value-of select="@compname" /></b></span>
      </h2>
      <div>
	<ul>
	  <xsl:apply-templates select="rest_file_extension"/>
	  <xsl:apply-templates select="hist_file_extension"/>
	  <xsl:apply-templates select="rest_history_varname"/>
	  <xsl:apply-templates select="rpointer"/>
	</ul>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="rest_file_extension">
    <li>Restart File Extension &lt;rest_file_extension&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="hist_file_extension">
    <li>History File Extension &lt;hist_file_extension&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="rest_history_varname">
    <li>Restart History Variable Name &lt;rest_history_varname&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="rpointer">
    rpointer File Specifications &lt;rpointer&gt; :
    <ul>
      <xsl:apply-templates select="rpointer_file"/>
      <xsl:apply-templates select="rpointer_content"/>
    </ul>
  </xsl:template>

  <xsl:template match="rpointer_file">
    <li>rpointer Filename &lt;rpointer_file&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="rpointer_content">
    <li>rpointer File Content &lt;rpointer_content&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

</xsl:stylesheet>
