<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="history_fields">
  <head>
    <title>CLM History Fields</title>
  </head>
  <body>
<hr/>
    <h2>Definition of CLM history variables</h2>
    <p>Included in the table are the following pieces of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Long name description.</li>
    <li>units</li>
    </ul>

    <table border="1" cellpadding="10">
    <caption>CLM History Fields</caption>
      <tr>
      <th>Name</th>
      <th>Long-name</th>
      <th>Units</th>
      </tr>
      <xsl:for-each select="*">
      <xsl:sort select="@name"/>
         <tr>
           <td><font color="#ff0000"><xsl:value-of select="@name"/></font></td>
           <td><xsl:value-of select="@long_name"/></td>
           <td><xsl:value-of select="@units"/></td>
         </tr>
      </xsl:for-each>
    </table>
<hr/>

  </body>
</xsl:template>

</xsl:stylesheet>
