<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="xml" version="1.0" encoding="UTF-8" indent="yes"/>

<xsl:strip-space elements="*"/>

<xsl:template match="/">
<!-- create the same node hierarchy as the main xml file -->
<scenario>
	<world>
		<!-- do for each row in the imported xml file or each row of the text file -->
		<xsl:for-each select="Import/Row">
			<!-- add region name as an attribute -->
			<region> <xsl:attribute name="name"><xsl:value-of select="Region"/></xsl:attribute>
				<!-- add sector name as an attribute -->
				<demandsector><xsl:attribute name="name"><xsl:value-of select="Sector"/></xsl:attribute>
				  	<!-- aeei is the new variable to be merged into the main xml file -->
					<aeei><xsl:attribute name="year">1975</xsl:attribute>
						<xsl:value-of select="y1975"/></aeei>
					<aeei><xsl:attribute name="year">1990</xsl:attribute>
						<xsl:value-of select="y1990"/></aeei>
					<aeei><xsl:attribute name="year">2005</xsl:attribute>
						<xsl:value-of select="y2005"/></aeei>
					<aeei><xsl:attribute name="year">2020</xsl:attribute>
						<xsl:value-of select="y2020"/></aeei>
					<aeei><xsl:attribute name="year">2035</xsl:attribute>
						<xsl:value-of select="y2035"/></aeei>
					<aeei><xsl:attribute name="year">2050</xsl:attribute>
						<xsl:value-of select="y2050"/></aeei>
					<aeei><xsl:attribute name="year">2065</xsl:attribute>
						<xsl:value-of select="y2065"/></aeei>
					<aeei><xsl:attribute name="year">2080</xsl:attribute>
						<xsl:value-of select="y2080"/></aeei>
					<aeei><xsl:attribute name="year">2095</xsl:attribute>
						<xsl:value-of select="y2095"/></aeei>
				</demandsector>				
			</region>
		</xsl:for-each>
	</world>
</scenario>
</xsl:template>
</xsl:stylesheet>
