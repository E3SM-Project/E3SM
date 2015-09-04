<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl  ="http://www.w3.org/1999/XSL/Transform"   version="2.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_compset">
  <head>
    <title>Configuration Component Sets</title>
  </head>
  <body>
    <h2>Configuration Component Sets</h2>

    <table BORDER="1" CELLPADDING="10">
      <th>ShortName (alias)</th>
      <th>Long name Description</th>
      <tr>
      <th colspan="3">A (All Data Models)</th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'A')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PA')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RA')]"/>
      <tr>
      <th colspan="3">B (All Active Models)</th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'B')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PB')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RB')]"/>
      <tr>
      <th colspan="3">C (Standalone POP)</th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'C')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PC')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RC')]"/>
      <tr>
      <th colspan="3">D (Active sea-ice and ocean, data atmosphere and stub land)</th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'D')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PD')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RD')]"/>
      <tr>
      <th colspan="3">E (Active land and atmosphere with slab ocean) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'E')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PE')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RE')]"/>
      <tr>
      <th colspan="3">F (Active land and atmosphere with data ocean) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'F')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PF')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RF')]"/>
      <tr>
      <th colspan="3">G (Active sea-ice and ocean, data atmosphere and land) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'G')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PG')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RG')]"/>
      <tr>
      <th colspan="3">H (Active sea-ice and ocean, data atmosphere and stub land) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'H')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PH')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RH')]"/>
      <tr>
      <th colspan="3">I (Standalone CLM) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'I')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PI')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RI')]"/>
      <tr>
      <th colspan="3">S (All stub models with dead atmosphere model) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'S')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PS')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RS')]"/>
      <tr>
      <th colspan="3">T (All stub models with data land model) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'T')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PT')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RT')]"/>
      <tr>
      <th colspan="3">U (stub models with data land model and active runoff [RTM]) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'U')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PU')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RU')]"/>
      <tr>
      <th colspan="3">X (All dead models) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'X')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'PX')]"/>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'RX')]"/>
      <tr>
      <th colspan="3">GLC (All compsets from above that have the active land-ice component) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[contains(@sname,'_GLC')]"/>
      <tr>
      <th colspan="3">P (All compsets using WRF regional atmosphere model with CLM) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'P')]"/>
      <tr>
      <th colspan="3">R (All compsets using WRF regional atmosphere model with VIC 
                      regional land model) </th>
      </tr>
      <xsl:apply-templates select="COMPSET[starts-with(@sname,'R')]"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="COMPSET">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@sname"/></font>
    (<xsl:value-of select="@alias"/>)
    </td>
    <td>
        <p><xsl:value-of select="text()"/></p>
    </td>
  </tr>
</xsl:template>


</xsl:stylesheet>
