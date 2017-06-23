<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" indent="yes" 
    doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"
    doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN" />

  <xsl:template match="//file">
    <html>
      <head>
	<title>CESM CASEROOT env_mach_specific.xml file</title>
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
	  <h1>CESM CASEROOT env_mach_specific.xml file</h1>
	  <p style="font-size: 0.9em;">
	    Model Version: CESM2.0<br/>
	    HTML created on: 2017-06-21
	  </p>

	  <h3>Description</h3>
	  <p>
	    The CASEROOT scripts read the XML settings in the env_mach_specific.xml file <br/>
	    in order to determine the machine specific system settings including necessary<br/>
	    dependency modules, paths, libraries, environment settings and command arguments.
	  </p>

	  <h3>Interface Tools</h3>
	  <p>
	    CESM2.0 supports 
	    the <a href="http://www.cesm.ucar.edu/models/cesm2.0/cesm/xmlquery.html">xmlquery</a> 
	    and <a href="http://www.cesm.ucar.edu/models/cesm2.0/cesm/xmlchange.html">xmlchange</a> script tools
	    to modify some of the settings in this file. <br/>
	    Changes to these settings may be made prior to calling <b>case.setup</b>.<br/>
	    Subsequent changes require <b>case.setup --reset</b> followed by <b>case.build --clean</b>
	    and <b>case.build</b>.
	  </p>
	  <h4>Examples</h4>

	  <h3>Navigation Instuctions</h3>
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
      <h2><span id="variable">File Description Header &lt;header&gt;</span></h2>
      <div>
	<pre>
	  <xsl:value-of select="." />
	</pre>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="module_system">
    <div class="accordionItem">
      <h2>Module System  <span id="variable"><b> &lt;module_system type="<xsl:value-of select="@type" />"&gt; : </b></span></h2>
      <div>
	  <xsl:apply-templates select="init_path"/>
	  <xsl:apply-templates select="cmd_path"/>
	  <xsl:apply-templates select="modules"/>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="init_path">
    <li>Initialization Path for Script Language &lt;init_path lang="<xsl:value-of select="@lang" />"&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="cmd_path">
    <li>Command Path for Script Language &lt;cmd_path lang="<xsl:value-of select="@lang" />"&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="modules">
    <dl>
      <dt>System Modules &lt;modules&gt; :</dt>
      <dd>
	<ul>
	  <li>Attributes :
	    <ul>
	      <li>compiler = <b><xsl:value-of select="@compiler" /></b></li>
	      <li>debug = <b><xsl:value-of select="@debug" /></b></li>
	      <li>mpilib = <b><xsl:value-of select="@mpilib" /></b></li>
	    </ul>
	  </li>
	  <xsl:apply-templates select="command"/>
	</ul>
      </dd>
    </dl>
  </xsl:template>

  <xsl:template match="command">
    <li>Module Command &lt;command name=<b><xsl:value-of select="@name" /></b>&gt; <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="environment_variables">
    <div class="accordionItem">
      <h2>Environment Variables <span id="variable">&lt;environment_variables&gt; :</span></h2>
      <div>
	<ul>
	  <xsl:apply-templates select="env"/>
	</ul>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="env">
    <li>&lt;env name=<b><xsl:value-of select="@name" /></b>&gt; <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="mpirun">
    <div class="accordionItem">
      <h2>MPI Run Commands <span id="variable">&lt;mpirun&gt; :</span></h2>
      <div>
	<ul>
	  <li>Attributes :
	    <ul>
	      <li>mpilib = <b><xsl:value-of select="@mpilib" /></b></li>
	      <li>unit_testing = <b><xsl:value-of select="@unit_testing" /></b></li>
	    </ul>
	  </li>
	  <xsl:apply-templates select="arguments"/>
	  <xsl:apply-templates select="executable"/>
	</ul>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="arguments">
    <li>MPI Run Command Arguments &lt;arguments&gt;
      <ul>
	  <xsl:apply-templates select="arg"/>
      </ul>
    </li>
  </xsl:template>

  <xsl:template match="arg">
    <li>&lt;arg name=<b><xsl:value-of select="@name" /></b>&gt; <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="executable">
    <li>MPI Run Executable Path &lt;executable&gt; <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="entry">
    <div class="accordionItem">
      <h2>Entry ID/Value Definitions <span id="variable">&lt;entry&gt; :</span></h2>
      <div>
	<ul>
	  <li>Attributes :
	    <ul>
	      <li>id = <b><xsl:value-of select="@id" /></b></li>
	      <li>value = <b><xsl:value-of select="@value" /></b></li>
	    </ul>
	  </li>
	  <xsl:apply-templates select="type"/>
	  <xsl:apply-templates select="values"/>
	  <xsl:apply-templates select="valid_values"/>
	  <xsl:apply-templates select="desc"/>
	</ul>
      </div>
    </div>
  </xsl:template>

  <xsl:template match="type">
    <li>Data Type &lt;type&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="values">
    <dl>
      <dt>Values &lt;values&gt; :</dt>
      <dd>
	<ul>
	  <xsl:apply-templates select="value"/>
	</ul>
      </dd>
    </dl>
  </xsl:template>

  <xsl:template match="value">
    <li>Value &lt;value component=<b><xsl:value-of select="@component" /></b>&gt; <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="valid_values">
    <li>Valid Values &lt;valid_values&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

  <xsl:template match="desc">
    <li>Description &lt;desc&gt; : <b><xsl:value-of select="." /></b></li>
  </xsl:template>

</xsl:stylesheet>
