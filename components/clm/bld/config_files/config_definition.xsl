<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_definition">
  <head>
    <title>CLM Configuration Definition</title>
  </head>
  <body>
    <h2>CLM Configuration Definition</h2>

    <table border="1" cellpadding="10">
    <caption><font size="larger"><bold>CLM Physics Configurations</bold></font></caption>
    <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Value</th>
      <th>Description</th>
    </tr>
    <tr>
      <th>Valid Values</th>
    </tr>
      <xsl:apply-templates select="entry[@category='physics']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption><font size="larger"><bold>CLM Biogeochemistry Configurations</bold></font></caption>
    <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Value</th>
      <th>Description</th>
    </tr>
    <tr>
      <th>Valid Value</th>
    </tr>
      <xsl:apply-templates select="entry[@category='bgc']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption><font size="larger"><bold>Configuration Directories</bold></font></caption>
    <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Value</th>
      <th>Description</th>
    </tr>
    <tr>
      <th>Valid Value</th>
    </tr>
      <xsl:apply-templates select="entry[@category='directories']"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td rowspan="2"><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td rowspan="2"><xsl:value-of select="@value"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
  <tr>
    <td><b>Valid values: </b> <xsl:value-of select="@valid_values"/></td>
  </tr>
</xsl:template>


</xsl:stylesheet>
