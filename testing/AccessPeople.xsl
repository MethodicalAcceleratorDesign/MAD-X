<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:output omit-xml-declaration="yes"/>
	<xsl:output method="text"/>

<xsl:param name="what" >Undefined</xsl:param><!-- first argument supplied on the command-line when invoking the xslt processor-->
<xsl:param name="who">Undefined</xsl:param><!-- second (optional) argument supplied on the command-line when invoking the xslt processor-->

<xsl:template match="/">
<xsl:choose>

	<xsl:when test="$what='Undefined'">Error, this stylesheet expects _what_ argument
	</xsl:when>

	<xsl:when test="$what='email'">
		<xsl:choose>
			<xsl:when test="$who='Undefined'">Error, this stylesheet expects _who_ argument
			</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="//person[@login=$who]/@e-mail"/>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:when>
	
</xsl:choose>
</xsl:template>
</xsl:stylesheet>
