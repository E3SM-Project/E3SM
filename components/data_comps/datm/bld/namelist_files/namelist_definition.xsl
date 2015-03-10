<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="namelist_definition">
  <head>
    <title>CESM DATM (Data Atmosphere Model) Namelist Definition</title>
  </head>
  <body>
<p>
</p>
<hr/>
<p>
</p>
    <h2>DATM namelist variables</h2>
    <p>Note, these all would go into the datm.buildnml file:</p>
    <p>Included in the table are the following pieces of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Variable type (<code>char</code>, <code>integer</code>,
    <code>real</code>, or <code>logical</code>).  The type
    <code>char</code> has the length appended
    following an asterisk, e.g., <code>char*256</code>.  Variables that are
    arrays have their dimension specifier appended inside parentheses.  For
    example <code>char*1(6)</code> denotes a array of six
    <code>char*1</code> values.
    </li>
    <li>Variable description (includes information on defaults).</li>
    <li>Valid values (if restricted).</li>
    </ul>

    <table border="1" cellpadding="10">
    <caption>Data Atmosphere Model Streams File Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='streams']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>Data Atmosphere Model Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm']"/>
    </table>
    <p>
    <b>NOTE:</b> The following settings are NOT namelist settings, but
     settings used by datm-build-namelist to figure out actual namelist
     settings.
    </p>
    <table border="1" cellpadding="10">
    <caption>Data Atmosphere Model Internal Namelist Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datm_setting']"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td rowspan="2"><font color="#ff0000"><xsl:value-of select="@id"/></font></td>
    <td rowspan="2"><xsl:value-of select="@type"/></td>
    <td><xsl:apply-templates/></td>
  </tr>
  <tr>
    <td colspan="1"><xsl:if test="string-length(@valid_values)>0"><b>Valid Values: </b>
         <xsl:value-of select="@valid_values"/></xsl:if></td>
  </tr>
</xsl:template>

</xsl:stylesheet>
