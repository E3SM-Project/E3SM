<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="namelist_definition">
  <head>
    <title>CLM Namelist Definition</title>
  </head>
  <body>
<p>
</p>
<hr/>
<p>
</p>
    <h1>Definition of CLM namelist variables</h1>
    <p>We list all of the relevant namelist variables for CLM I cases. This includes
    CLM Namelist items as well as CLM build-namelist settings and namelist settings
    for CLM offline tools.</p>
<hr/>
    <h2>Definition of CLM namelist variables</h2>
    <p>Note, these all would go into the user_nl_clm file
    before configure):</p>
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
    <caption>CLM Namelist Physics Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_physics']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Lake Model Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_lake']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Permafrost Model Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_permafrost']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Methane Model Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_methane']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Vertical CN Model Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_vertcn']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Carbon Isotope Model Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_isotope']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Datasets</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description
      </th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='datasets']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist History output settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='history']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Restart settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_restart']"/>
    </table>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Performance Tuning</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='clm_performance']"/>
    </table>

<p>
</p>
<hr/>
<p>
</p>
    <h2>Command Line Options to CLM Build-namelist</h2>
    <p>Variables that are entered as options to build-namelist (but NOT used by 
        namelists in code). Most of these are options that could be added to
        CLM_BLDNML_OPTS. Included in the table are the following pieces 
        of information:</p>
    <ul>
    <li>Variable name.</li>
    <li>Type.</li>
    <li>Valid values.</li>
    <li>Variable description.</li>
    </ul>

    <table border="1" cellpadding="10">
    <caption>CLM Namelist Default Settings</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='default_settings']"/>
    </table>
<p>
</p>
<hr/>
<p>
</p>
    <h2>Namelist items for CLM Tools</h2>
    <p>These are namelist items that appear in the CLM Tools under models/lnd/clm/tools.
    </p>
    <table border="1" cellpadding="10">
    <caption>CLM mksurfdata</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='mksurfdata']"/>
    </table>
    <table border="1" cellpadding="10">
    <caption>CLM mkgriddata</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='mkgriddata']"/>
    </table>
    <table border="1" cellpadding="10">
    <caption>Miscellaneous CLM tools</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values</th>
      </tr>
      <xsl:apply-templates select="entry[@category='tools']"/>
    </table>

<hr/>

<h2>Namelist items for Driver Dry Deposition</h2>
    <table border="1" cellpadding="10">
    <caption>Driver Dry-Deposition Namelist Options</caption>
      <tr>
      <th rowspan="2">Name</th>
      <th rowspan="2">Type</th>
      <th>Description</th>
      </tr>
      <tr>
      <th colspan="1">Valid values, if restricted at all</th>
      </tr>
      <xsl:apply-templates select="entry[@category='dry-deposition']"/>
    </table>

<hr/>

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
