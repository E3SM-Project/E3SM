<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
  <head>
    <title>CESM DATM (Data Atmosphere Model) Namelist Defaults</title>
  </head>
  <body>
    <h2>Default Values for Namelist Variables</h2>
    <p>Included in the table are the following pieces of information:</p>
    <h3>Table headers include:</h3>
    <ul>
       <li>Name of variable</li>
       <li>Stream</li>
       <li>Grid</li>
    </ul>
    <h3>Miscellaneous items include:</h3>
    <ol>
       <li>DATM_MODE</li>
       <li>presaero_flag</li>
       <li>presaero_mode</li>
    </ol>

    <table border="1" cellpadding="10">
    <caption>Namelist Defaults</caption>
    <tr>
      <th rowspan="2">Name</th>
      <th>Stream</th>
      <th>Horz. Grid</th>
      <th>Miscellaneous</th>
    </tr>
    <tr>
      <th colspan="5">Default Value for this Configuration</th>
    </tr>
      <xsl:for-each select="namelist_defaults/*">
      <xsl:sort select="name()"/>
      <tr>
        <td rowspan="2"><font color="#ff0000">
        <xsl:value-of select="name()"/>
        </font></td>
        <td>
        <xsl:choose>
        <xsl:when test="string-length(@stream)>0">
             <xsl:value-of select="@stream"/>
        </xsl:when>
        <xsl:otherwise>
             All streams
        </xsl:otherwise>
        </xsl:choose>
        </td>
        <td>
        <xsl:choose>
        <xsl:when test="string-length(@grid)>0">
             <xsl:value-of select="@grid"/>
        </xsl:when>
        <xsl:otherwise>
             All res
        </xsl:otherwise>
        </xsl:choose>
        </td>
        <td>
        <xsl:if test="string-length(@datm_mode)>0">
        DATM_MODE=<xsl:value-of select="@datm_mode"/>
        </xsl:if>
        <xsl:if test="string-length(@presaero_flag)>0">
        presaero_flag=<xsl:value-of select="@presaero_flag"/>
        </xsl:if>
        <xsl:if test="string-length(@presaero_mode)>0">
        presaero_mode=<xsl:value-of select="@presaero_mode"/>
        </xsl:if>
        </td>
      </tr>
      <tr>
        <td colspan="5"><b>Value: </b><xsl:value-of select="."      /></td>
      </tr>
      </xsl:for-each>
    </table>

  </body>

  </html>
</xsl:template>

</xsl:stylesheet>
