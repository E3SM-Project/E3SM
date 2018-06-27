<?xml version='1.0' ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xmlns="http://www.w3.org/2000/xmlns/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
	<xsl:template match="/">
		<scenario>
			<xsl:value-of select="scenario/@date"/>
			<xsl:value-of select="scenario/@name"/>
			<xsl:value-of select="scenario/@xmlns:xsi"/>
			<xsl:value-of select="scenario/@xsi:noNamespaceSchemaLocation"/>
			<xsl:value-of select="scenario/summary"/>
			<xsl:value-of select="scenario/modeltime"/>
			<xsl:value-of select="scenario/world"/>
		</scenario>
	</xsl:template>
</xsl:stylesheet><!-- Stylus Studio meta-information - (c) 2004-2005. Progress Software Corporation. All rights reserved.
<metaInformation>
<scenarios ><scenario default="yes" name="base_input.xml" userelativepaths="yes" externalpreview="no" url="base_input.xml" htmlbaseurl="" outputurl="" processortype="internal" useresolver="yes" profilemode="0" profiledepth="" profilelength="" urlprofilexml="" commandline="" additionalpath="" additionalclasspath="" postprocessortype="none" postprocesscommandline="" postprocessadditionalpath="" postprocessgeneratedext="" validateoutput="no" validator="internal" customvalidator=""/></scenarios><MapperMetaTag><MapperInfo srcSchemaPathIsRelative="yes" srcSchemaInterpretAsXML="no" destSchemaPath="" destSchemaRoot="" destSchemaPathIsRelative="yes" destSchemaInterpretAsXML="no" ><SourceSchema srcSchemaPath="file:///c:/work_area/Code/objects/workspace&#x2D;vintage2/cvs/objects/model_data/base/base_input.xml" srcSchemaRoot="scenario" AssociatedInstance="" loaderFunction="document" loaderFunctionUsesURI="no"/></MapperInfo><MapperBlockPosition><template match="/"></template></MapperBlockPosition><TemplateContext></TemplateContext><MapperFilter side="source"></MapperFilter></MapperMetaTag>
</metaInformation>
-->