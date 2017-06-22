<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" indent="yes" 
    doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"
    doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN" />

  <xsl:template match="//file">
    <html>
      <head>
	<title>CESM CASEROOT env_batch.xml file</title>
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
	  <h1>CESM CASEROOT example env_batch.xml file</h1>
	  <p style="font-size: 0.9em;">
	    Model Version: CESM2.0<br/>
	    HTML created on: 2017-06-21
	  </p>
  	  <p>
	    This page contains an example env_batch.xml file read by the case.submit script.
	  </p>
	  <p>
	    The case.submit reads these XML settings in order to<br/>
	    determine the machine specific batch job submission arguments to be used for<br/>
	    job group. CESM2.0 supports 3 job groups; <i>case.run</i>, <i>case.test</i>, and <i>case.st_archive</i>.<br/>
	    The <a href="http://www.cesm.ucar.edu/models/cesm2.0/cesm/xmlquery.html">xmlquery</a> 
	    and <a href="http://www.cesm.ucar.edu/models/cesm2.0/cesm/xmlchange.html">xmlchange</a> script tools<br/>
	    accept the <b>--subgroup</b> command line argument to define one of the 3 job groups.
	  </p>
	  <p>
	    Varible names delimited by <pre>{{ }}</pre> are substituted at run-time and should not be modified.
	  </p>
	  <p>
	  Clicking on the blue link will display additional descriptive information.<br/>  
	  Click on the "Show Details" button and then cntl+F key to search for specific strings in this file.
	  </p>

	  <button class="detailsBtn" onclick="showAll();">Show Details</button>
	  <button class="detailsBtn" onclick="hideAll();">Hide Details</button>

	  <br/><br/>
	  <xsl:apply-templates />
	</div>
      </body>
    </html>
  </xsl:template>

  <xsl:template match="header">
    <div class="accordionItem">
      <h2><span id="variable">Header &lt;header&gt;</span></h2>
      <div>
	<pre>
	  <xsl:value-of select="//file/header" />
	</pre>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="group[@id='config_batch']">
    <div class="accordionItem">
      <h2>Default Batch Group Settings <span id="variable"><b>&lt;group id="config_batch"&gt; : </b></span></h2>
      <div>
	  <xsl:apply-templates select="entry"/>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="batch_system[@type='pbs'][1]">
    <div class="accordionItem">
      <h2>Default Batch System Settings <span id="variable"><b>&lt;batch_system type="pbs"&gt; : </b></span></h2>
      <div>
	<ul>
	  <xsl:apply-templates select="batch_query"/>
	  <xsl:apply-templates select="batch_submit"/>
	  <xsl:apply-templates select="batch_directive"/>
	  <xsl:apply-templates select="jobid_pattern"/>
	  <xsl:apply-templates select="depend_string"/>
	  <xsl:apply-templates select="walltime_format"/>
	  <xsl:apply-templates select="submit_args"/>
	  <xsl:apply-templates select="directives"/>
	</ul>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="batch_system[@MACH='cheyenne' and @type='pbs']">
    <div class="accordionItem">
      <h2>Machine Specific Batch System Settings <span id="variable"><b>&lt;batch_system MACH="cheyenne" type="pbs"&gt; : </b></span></h2>
      <div>
	<ul>
	  <xsl:apply-templates select="directives"/>
	  <xsl:apply-templates select="queues"/>
	</ul>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="group[@id='case.run']">
    <div class="accordionItem">
      <h2>Group <i>case.run</i> Settings <b><span id="variable">&lt;group id="case.run"&gt; : </span></b></h2>
      <div>
	  <xsl:apply-templates select="entry"/>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="group[@id='case.test']">
    <div class="accordionItem">
      <h2>Group <i>case.run</i> Settings <b><span id="variable">&lt;group id="case.test"&gt; : </span></b></h2>
      <div>
	  <xsl:apply-templates select="entry"/>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="group[@id='case.st_archive']">
    <div class="accordionItem">
      <h2>Group <i>case.st_archive</i> Settings <b><span id="variable">&lt;group id="case.st_archive"&gt; : </span></b></h2>
      <div>
	  <xsl:apply-templates select="entry"/>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="entry">
    <dl>
      <dt>entry &lt;entry&gt; :</dt>
      <dd>
	<ul>
	  <li> attributes :
	    <ul>
	      <li>id = <b><xsl:value-of select="@id" /></b></li>
	      <li>value = <b><xsl:value-of select="@value" /></b></li>
	    </ul>
	  </li>
	  <xsl:apply-templates select="type"/>
	  <xsl:apply-templates select="valid_values"/>
	  <xsl:apply-templates select="desc"/>
	</ul>
      </dd>
    </dl>
  </xsl:template>

  <xsl:template match="type">
    <li>Data Type &lt;type&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="valid_values">
    <li>Valid Values &lt;valid_values&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="desc">
    <li>Description &lt;desc&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="batch_query">
    <li>Batch Query Command &lt;batch_query&gt; : <b><xsl:value-of select="." /></b>
      <ul>
	<li> attributes :
	  <ul>
	    <li>args = <b><xsl:value-of select="@args" /></b></li>
	  </ul>
	</li>
      </ul>
    </li>
  </xsl:template>

  <xsl:template match="batch_submit">
    <li>Batch Submission Command &lt;batch_submit&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="batch_directive">
    <li>Batch Directive Prefix &lt;batch_directive&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="jobid_pattern">
    <li>Job ID Regular Expression Pattern &lt;jobid_pattern&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="depend_string">
    <li>Job Dependency Command &lt;depend_string&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="walltime_format">
    <li>Job Walltime Format &lt;walltime_format&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="submit_args">
    Batch Job Submission Arguments &lt;submit_args&gt; :
    <ul>
      <xsl:apply-templates select="arg"/>
    </ul>
  </xsl:template>

  <xsl:template match="arg">
    <li>Batch Argument &lt;arg&gt; : 
      <ul>
	<li>Attributes :
	  <ul>
	    <li>flag = <b><xsl:value-of select="@flag" /></b></li>
	    <li>name = <b><xsl:value-of select="@name" /></b></li>
	  </ul>
	</li>
      </ul>
    </li>
  </xsl:template>

  <xsl:template match="directives">
    Batch Job Submission Directives &lt;direcitves&gt; :
    <ul>
      <xsl:apply-templates select="directive"/>
    </ul>
  </xsl:template>

  <xsl:template match="directive">
    <li>&lt;directive&gt; : <b><xsl:value-of select="." /></b><br/>
      default=<b><xsl:value-of select="@default" /></b><br/>
    </li>
  </xsl:template>

  <xsl:template match="queues">
    Batch Submission Default Queue Settings &lt;queues&gt; : <br/>
    Note: each job group defines their own job submission settings.
    <ul>
      <xsl:apply-templates select="queue"/>
    </ul>
  </xsl:template>

  <xsl:template match="queue">
    <li>Queue Name &lt;queue&gt; : <b><xsl:value-of select="." /></b>
      <ul>
	<li>Attributes :
	  <ul>
	    <li>default = <b><xsl:value-of select="@default" /></b></li>
	    <li>jobmax = <b><xsl:value-of select="@jobmax" /></b></li>
	    <li>jobmin = <b><xsl:value-of select="@jobmin" /></b></li>
	    <li>walltimemax = <b><xsl:value-of select="@jobmin" /></b></li>
	  </ul>
	</li>
      </ul>
    </li>
  </xsl:template>

</xsl:stylesheet>
