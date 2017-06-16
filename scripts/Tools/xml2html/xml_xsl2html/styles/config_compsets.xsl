<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns="http://www.w3.org/1999/xhtml" version="1.0">
  <xsl:output encoding="UTF-8" indent="yes" method="xml" standalone="no" omit-xml-declaration="no"/>

  <xsl:template match="/">
    <html>
      <head>
	<title>CIME CESM Allactive config_compsets.xml</title>
	<!-- link type="text/css" href="/styles/cesm.css" rel="stylesheet"/ -->
      </head>

      <body>

	<h1>CIME CESM Allactive config_compsets.xml</h1>
	<h2>How to read/parse the component set long names.</h2>
	<pre>
	  <xsl:value-of select="compsets/help" />
	</pre>
	<hr/>

	<table bgcolor="LightGray" border="1">
	  <tbody>
	    <tr>
	      <th>Alias</th>
	      <th>Long Name</th>
	      <th>Comment</th>
	    </tr>
	    <xsl:for-each select="compsets/compset">
	    <tr>
	      <td>
		<xsl:value-of select="alias"/>
	      </td>
	      <td>
		<xsl:value-of select="lname"/>
	      </td>
	      <td>
		<xsl:if test="preceding-sibling::comment()[1]">
		  <xsl:value-of select="preceding-sibling::comment()[1]" />
		</xsl:if>
	      </td>
	    </tr>
	    </xsl:for-each>
	  </tbody>
	</table>
	<hr/>
	<xsl:apply-templates />
      </body>
    </html>
  </xsl:template>

  <hr/>
  <h2>Run specific settings based on matching component set and grid long names.</h2>

  <xsl:template match="compsets/entries/entry">
    <p>
      XML variable: <xsl:value-of select="@id" />
    </p>
    <xsl:apply-templates />
  </xsl:template>

  <xsl:template match="compsets/help" />

  <xsl:template match="compsets/compset" />

  <xsl:template match="//compsets/entries/entry/values/value">
    <ul>
      <li>Grid: <xsl:value-of select="@grid" /></li>
      <li>Compset: <xsl:value-of select="@compset" /></li>
      <li>Value: <xsl:value-of select="." /></li>
    </ul>
  </xsl:template>

</xsl:stylesheet>
