<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_definition">
  <head>
    <title>Configuration Definition</title>
  </head>
  <body>
    <h2>Configuration Definition</h2>

    <table BORDER="1" CELLPADDING="10">
      <th>Name</th>
      <th>Value</th>
      <th>Description</th>
      <xsl:apply-templates select="entry"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td><xsl:value-of select="@value"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
</xsl:template>


</xsl:stylesheet>
