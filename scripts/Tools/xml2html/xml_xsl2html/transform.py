import lxml.etree as ET

xml_filename = "./config_compsets.xml"
xsl_filename = "./styles/config_compsets.xsl"

dom = ET.parse(xml_filename)
xslt = ET.parse(xsl_filename)
transform = ET.XSLT(xslt)
newdom = transform(dom)
print(ET.tostring(newdom, pretty_print=True))
