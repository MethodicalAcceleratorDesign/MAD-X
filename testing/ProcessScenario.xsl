<?xml version="1.0"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:output omit-xml-declaration="yes"/>
	<xsl:output method="text"/>

<xsl:param name="what" >Undefined</xsl:param><!-- first argument supplied on the command-line when invoking the xslt processor-->
<xsl:param name="target">Undefined</xsl:param><!-- second (optional) argument supplied on the command-line when invoking the xslt processor-->

<xsl:template match="/">
<xsl:choose>

	<xsl:when test="$what='list_targets'">
		<xsl:for-each select="//target">
			<xsl:value-of select="./@name"/><xsl:text>
</xsl:text>
		</xsl:for-each>

	</xsl:when>
	<xsl:when test="$what='Undefined'">Error, this stylesheet expects _what_ argument
	</xsl:when>

	<xsl:when test="$what='list_tests'">
		<xsl:choose>
			<xsl:when test="$target='Undefined'">Error, this stylesheet expects _target_ argument
			</xsl:when>
			<xsl:otherwise>
<!-- in examples, the redirection is "> !" but the command failed to produce the output file when invoked from a perl-script. Solved the problem by removing the "!" and only retaining the ">" -->
				<xsl:for-each select="//target[@name=$target]/test"><xsl:choose><xsl:when test="./@program">./<xsl:value-of select="./@program"/></xsl:when><xsl:otherwise>./madx</xsl:otherwise></xsl:choose> &lt; <xsl:value-of select="./@input"/> &gt; <xsl:text> </xsl:text> <xsl:value-of select="./@output"/> <xsl:if test="./@subdir">, subdirectory=<xsl:value-of select="./@subdir"/></xsl:if><xsl:text>
</xsl:text>
				</xsl:for-each>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:when>
	

</xsl:choose>

</xsl:template>


</xsl:stylesheet>
