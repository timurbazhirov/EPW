<?xml version="1.0" encoding="ISO-8859-1"?>

<!--
    ***
    *** THIS FILE IS a XSL STYLESHEET FOR TRANSFORMING INPUT_*.xml to INPUT_*.html
    ***
-->

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">    

  <!--<xsl:strip-space elements="*"/>-->
  <xsl:output method="html"/>  

  <!-- *** ROOT *** -->

  <xsl:template match="/input_description">
    <html>
      <head>
	<xsl:comment> *** FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST *** </xsl:comment>
	<meta http-equiv="Content-Style-Type" CONTENT="text/css" />
      </head>
      <body style="background-color: white; font-family: arial, helvetica, sans-serif; font-size: 14px; padding-left: 10; width: 700; ">
	<a name="__top__"></a>
	<table style="border-width: 0 width: 100% text-align: left; vertical-align: top; background: #00395a;">
	  <tr>
	    <th style="margin: 3 3 3 10; background: #00395a; color: #ffffee; ">
	      <h1 style="margin: 5 10 10 15; text-align: left;"> Input File Description </h1>
	      <h2 style="margin: 5 10 10 15; text-align: left;"> Program: <xsl:value-of select="@program"/> / <xsl:value-of select="@package"/> / <xsl:value-of select="@distribution"/></h2>
	    </th>
	  </tr>
	  <tr><td style="padding: 10 3 3 3; background: #ffffff; color: #222222; ">
	    <xsl:apply-templates/>
	  </td></tr>
	</table>
	<blockquote>
	  <small>
	    This file has been created by helpdoc utility.
	  </small>
	</blockquote>
      </body>
    </html>
  </xsl:template>

  <!--  *** TOC ***  -->

  <xsl:template match="toc">
    <blockquote>
      <h3>TABLE OF CONTENTS</h3>
      <blockquote>	
	<xsl:for-each select="//intro">
	  <p><a href="#{generate-id(.)}">INTRODUCTION</a></p>
	</xsl:for-each>
	<xsl:for-each select="//namelist | //card | //linecard | //section">
	  <p>
	    <xsl:choose>
	      <xsl:when test="name(.)='linecard'">
		<a href="#{generate-id(.)}">Line-of-input:</a><xsl:text> </xsl:text>
		<xsl:for-each select=".//var | .//list">
		  <xsl:if test="name(.)='var' or name(.)='dimension'">
		    <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a> 
		    <xsl:if test="not(position()=last())">
		      <xsl:text> | </xsl:text>
		    </xsl:if>
		  </xsl:if>
		  <xsl:if test="name(.)='list'">
		    <a href="#{generate-id(.)}"><xsl:value-of select="format"/></a> 
		    <xsl:if test="not(position()=last())">
		      <xsl:text> | </xsl:text>
		    </xsl:if>
		  </xsl:if>
		</xsl:for-each>
	      </xsl:when>
	      <xsl:when test="name(.)='section'">
		<a href="#{generate-id(.)}"><xsl:value-of select="@title"/></a>		
		<xsl:apply-templates select="subsection" mode="toc"/>
	      </xsl:when>
	      <xsl:otherwise>
		<a href="#{generate-id(.)}">
		  <xsl:if test="name(.)='namelist'">&amp;</xsl:if><xsl:value-of select="@name"/>
		</a>
		<blockquote>
		  <xsl:for-each select=".//var | .//dimension | .//list | .//col | .//row">
		    <xsl:if test="name(.)='var' or name(.)='dimension'">
		      <xsl:if test="info != '' or 
				    status != '' or 
				    see    != '' or 
				    ../../vargroup/info != '' or 
				    ../../dimensiongroup/info != ''">
			<a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a> 
			<xsl:if test="not(position()=last())">
			  <xsl:text> | </xsl:text>
			</xsl:if>
		      </xsl:if>
		    </xsl:if>
		    <xsl:if test="name(.)='list'">
		      <xsl:if test="info != '' or status != '' or see != ''">
			<a href="#{generate-id(.)}"><xsl:value-of select="format"/></a> 
			<xsl:if test="not(position()=last())">
			  <xsl:text> | </xsl:text>
			</xsl:if>
		      </xsl:if>
		    </xsl:if>
		    <xsl:if test="name(.)='col'">
		      <xsl:if test="info != '' or status != '' or see != ''">
			<a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
			<xsl:if test="not(position()=last())">
			  <xsl:text> | </xsl:text>
			</xsl:if>
		      </xsl:if>
		      <xsl:if test="ancestor::colgroup/info != ''">
			<a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
			<xsl:if test="not(position()=last())">
			  <xsl:text> | </xsl:text>
			</xsl:if>
		      </xsl:if>
		    </xsl:if>
		    <xsl:if test="name(.)='row'">
		      <xsl:if test="info != '' or status != '' or see != ''">
			<a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
			<xsl:if test="not(position()=last())">
			  <xsl:text> | </xsl:text>
			</xsl:if>
		      </xsl:if>
		      <xsl:if test="ancestor::rowgroup/info != ''">
			<a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
			<xsl:if test="not(position()=last())">
			  <xsl:text> | </xsl:text>
			</xsl:if>
		      </xsl:if>
		    </xsl:if>
		  </xsl:for-each>
		</blockquote>
	      </xsl:otherwise>
	    </xsl:choose>
	  </p>
	</xsl:for-each>
      </blockquote>
    </blockquote>
  </xsl:template>

  <xsl:template match="subsection" mode="toc">
    <blockquote>
      <a href="#{generate-id(.)}"><xsl:value-of select="@title"/></a>
      <xsl:apply-templates select="subsubsection" mode="toc"/>
    </blockquote>
  </xsl:template>
  <xsl:template match="subsubsection" mode="toc">
    <blockquote>
      <a href="#{generate-id(.)}"><xsl:value-of select="@title"/></a>
    </blockquote>
  </xsl:template>

  <!-- back-to-top -->

  <xsl:template name="back_to_top">
    <div align="right">[<a href="#__top__">Back to Top</a>]</div>
  </xsl:template>

  <!--    *** INTRO ***  -->

  <xsl:template match="intro">
    <blockquote>
      <a name="{generate-id()}">
	<h3>INTRODUCTION</h3>
      </a>
      <blockquote>
	<pre>
	  <xsl:value-of select="."/>
	</pre>
      </blockquote>
    </blockquote>
  </xsl:template>


  <!--    *** GROUP ***  -->

  <xsl:template match="group">
    <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
      <tr><td>
	<xsl:apply-templates/>
      </td></tr>
    </table>
  </xsl:template>

  
  <!--    *** NAMELIST ***  -->

  <xsl:template match="namelist">
    <a name="{generate-id(.)}">
      <table border="0" width="100%" style="margin-bottom: 20; ">
	<tr>
	  <th bgcolor="#ddcba6">
	    <h2 style="margin: 10 10 10 15; text-align: left;"> Namelist: <xsl:value-of select="@name"/> </h2>
	  </th>
	</tr>
	<tr><td style="text-align: left; background: #ffebc6; padding: 5 5 5 30; ">	    
	  <table style="border-color: #505087; border-style: solid; border-width: 0; margin-bottom: 10; table-layout: auto; width: 700; ">
	    <tbody>
	      <tr><td>
		<xsl:apply-templates/>
	      </td></tr>
	    </tbody>
	  </table>
	</td></tr>
      </table>
    </a>      
  </xsl:template>


  <!--    *** CARD *** -->

  <xsl:template match="card">
    <a name="{generate-id(.)}">
      <table border="0" width="100%" style="margin-bottom: 20; ">
	<tr>
	  <th bgcolor="#ddcba6">
	    <h2 style="margin: 10 10 10 15; text-align: left;"> 
	    Card: <xsl:value-of select="@name"/>
	    <xsl:choose>
	      <xsl:when test="flag/@use = 'optional'">
		<xsl:text> { </xsl:text> 
		<xsl:value-of select="flag/enum" /> 
		<xsl:text> } </xsl:text>
	      </xsl:when>
	      <xsl:otherwise>
		<xsl:text> </xsl:text><xsl:value-of select="flag/enum" />
	      </xsl:otherwise>
	    </xsl:choose>
	    </h2>
	  </th>
	</tr>
	
	<tr><td style="text-align: left; background: #ffebc6; padding: 5 5 5 30; ">	    
	  <table style="border-color: #505087; border-style: solid; border-width: 0; margin-bottom: 10; table-layout: auto; width: 700; ">
	    <tbody>	      	      
	      <tr><td>
		<xsl:apply-templates select="syntax | if | choose | label | message"/>
	      </td></tr>
	      <tr><td>
		<h3>Description of items:</h3>
		<blockquote>
		  <pre><xsl:value-of select="flag/info"/></pre>
		  <xsl:apply-templates select="descendant::vargroup | descendant::var | descendant::list | descendant::table" mode="card_description"/>
		</blockquote>
	      </td></tr>
	    </tbody>
	  </table>
	</td></tr>
      </table>
    </a>      
  </xsl:template>

  <!-- card/syntax -->

  <xsl:template match="syntax">
    <h3>Syntax:</h3>
    <blockquote>
      <xsl:if test="boolean(ancestor::card/@nameless) = false()">
	<b>
	  <xsl:value-of select="ancestor::card/@name"/>
	  <xsl:choose>
	    <xsl:when test="normalize-space(@flag) = ''">
	      <xsl:message>empty -flag</xsl:message>
	      <xsl:choose>
		<xsl:when test="ancestor::card/flag/@use = 'optional'">
		  <xsl:text> { </xsl:text> 
		  <xsl:value-of select="ancestor::card/flag/enum" /> 
		  <xsl:text> } </xsl:text>
		</xsl:when>
		<xsl:when test="ancestor::card/flag/@use = 'conditional'">
		  <xsl:text> [ </xsl:text> 
		  <xsl:value-of select="ancestor::card/flag/enum" /> 
		  <xsl:text> ] </xsl:text>
		</xsl:when>
		<xsl:otherwise>
		  <xsl:text> </xsl:text><xsl:value-of select="ancestor::card/flag/enum" />
		</xsl:otherwise>
	      </xsl:choose>
	    </xsl:when>
	    <xsl:otherwise>
	      <xsl:message>non-empty -flag; <xsl:value-of select="name(.)"/>;<xsl:value-of select="@flag"/>;</xsl:message>
	      <xsl:text> </xsl:text><xsl:value-of select="@flag" /> 
	    </xsl:otherwise>
	  </xsl:choose>
	</b>
	<br/>
      </xsl:if>
      <xsl:apply-templates select="table | line | optional | conditional | list" mode="syntax"/>
    </blockquote>
  </xsl:template>

  <!-- card//syntax//line -->

  <xsl:template match="line" mode="syntax">
    <xsl:apply-templates select="optional | conditional | var | keyword | vargroup | list" mode="syntax"/>
    <br/>
  </xsl:template>	


  <!-- card//syntax//optional -->      

  <xsl:template match="optional" mode="syntax">
    <!--<div style="background: #eeeeee; color: #555555;">-->
    <xsl:text> { </xsl:text>
    <xsl:apply-templates select="line | var | list | keyword | table" mode="syntax"/>
    <xsl:text> } </xsl:text>	    
    <!--</div>-->
  </xsl:template>	      

  <!-- card//syntax//conditional -->      

  <xsl:template match="conditional" mode="syntax">
    <xsl:text> [ </xsl:text>
    <xsl:apply-templates select="line | var | list | keyword | table" mode="syntax"/>
    <xsl:text> ] </xsl:text>	    
  </xsl:template>	      

  <!-- card//syntax//keyword -->      

  <xsl:template match="keyword" mode="syntax">
    <b><xsl:value-of select="@name"/></b><xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//list -->

  <xsl:template match="list" mode="syntax">    
    <i>
      <xsl:choose>	
	<xsl:when test="info != ''">
	  <a href="#{generate-id(.)}"><xsl:value-of select="format"/></a>
	</xsl:when>
	<xsl:otherwise>
	  <xsl:value-of select="format"/>
	</xsl:otherwise>
      </xsl:choose>
    </i>
    <xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//var -->      

  <xsl:template match="var" mode="syntax">
    <xsl:message>var query = <xsl:value-of select="child::node()"/> </xsl:message>
    <i>
      <xsl:choose>	
	<xsl:when test="info != '' or status != '' or see != '' or ../../vargroup/info != ''">
	  <a href="#{generate-id(.)}"><xsl:value-of select="@name"/></a>
	</xsl:when>
	<xsl:otherwise>
	  <xsl:value-of select="@name"/>
	</xsl:otherwise>
      </xsl:choose>
    </i>
    <xsl:text>&#160;&#160;</xsl:text>
  </xsl:template>

  <!-- card//syntax//vargroup -->      

  <xsl:template match="vargroup" mode="syntax">
    <xsl:apply-templates select="var" mode="syntax"/>
  </xsl:template>

  <!-- card//syntax//table -->      

  <xsl:template match="table" mode="syntax">
    <a name="{generate-id(.)}">      
      <table>
	<xsl:apply-templates select="rows | cols" mode="syntaxTableMode"/>
      </table>
    </a>
  </xsl:template>

  <!-- sytntax//table/rows -->

  <xsl:template match="rows" mode="syntaxTableMode">
    <xsl:message>//card//syntax//table/rows</xsl:message>
    <tr>
      <xsl:call-template name="row">
	<xsl:with-param name="rowID"><xsl:value-of select="@start"/></xsl:with-param>
      </xsl:call-template>
    </tr> 
    <tr>
      <xsl:call-template name="row">
	<xsl:with-param name="rowID"><xsl:value-of select="number(@start+1)"/></xsl:with-param>
      </xsl:call-template>
    </tr> 
    <xsl:choose>
      <xsl:when test="number(@end) != @end">
	<tr><td colspan="2"><xsl:text>&#160;. . .</xsl:text></td></tr>
	<tr>
	  <xsl:call-template name="row">
	    <xsl:with-param name="rowID"><xsl:value-of select="@end"/></xsl:with-param>
	  </xsl:call-template>
	</tr> 
      </xsl:when>
      <xsl:otherwise>
	<xsl:if test="number(@end) > number(@start+2)">
	  <tr><td colspan="2"><xsl:text>&#160;. . .</xsl:text></td></tr>
	</xsl:if>
	<xsl:if test="number(@end) > number(@start+1)">
	  <tr>
	    <xsl:call-template name="row">
	      <xsl:with-param name="rowID"><xsl:value-of select="@end"/></xsl:with-param>
	    </xsl:call-template>
	  </tr>
	</xsl:if>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="row">
    <xsl:param name="rowID" select="1"/>   
    <xsl:message>//card//syntax//table//rows->rows(<xsl:value-of select="name(.)"/>)</xsl:message>
    <xsl:apply-templates select="col | colgroup | optional | conditional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="colgroup" mode="rowMode">
    <xsl:param name="rowID"/>
    <xsl:apply-templates select="col | optional | conditional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="optional" mode="rowMode">
    <xsl:param name="rowID"/>
    <td><xsl:text> { </xsl:text></td>
    <xsl:apply-templates select="col | colgroup | conditional | optional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
    <td><xsl:text> } </xsl:text></td>
  </xsl:template>

  <xsl:template match="conditional" mode="rowMode">
    <xsl:param name="rowID"/>
    <td><xsl:text> [ </xsl:text></td>
    <xsl:apply-templates select="col | colgroup | conditional | optional" mode="rowMode">
      <xsl:with-param name="rowID" select="$rowID"/>
    </xsl:apply-templates>
    <td><xsl:text> ] </xsl:text></td>
  </xsl:template>

  <xsl:template match="col" mode="rowMode">
    <xsl:param name="rowID"/>
    <td>  
      <xsl:text>&#160;</xsl:text>
      <i>
	<xsl:message>col query = <xsl:value-of select="child::node()"/> </xsl:message>

      <xsl:choose>	
	<xsl:when test="info != '' or status != '' or see != '' or ../../colgroup/info != ''">
	  <a href="#{generate-id(.)}"><xsl:value-of select="@name"/>(<xsl:value-of select="$rowID"/>)</a>
	</xsl:when>
	<xsl:otherwise>
	  <xsl:value-of select="@name"/>(<xsl:value-of select="$rowID"/>)
	</xsl:otherwise>
      </xsl:choose>
      </i>
      <xsl:text>&#160;</xsl:text>
    </td>  
  </xsl:template>
  
  <!-- syntax//table/cols -->

  <xsl:template match="cols" mode="syntaxTableMode">
    <xsl:message>//card//syntax//table/cols</xsl:message>
    <xsl:apply-templates select="row | rowgroup | optional | conditional" mode="colsMode">
      <xsl:with-param name="colsOptional"  select="false()"/>
      <xsl:with-param name="colsConditional" select="false()"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="row" mode="colsMode">
    <xsl:param name="colsOptional" select="false()"/>
    <xsl:param name="colsConditional" select="false()"/>
    <tr>
      <td align="right">
	<xsl:if test="$colsOptional    = true()"><xsl:text>{ &#160;</xsl:text></xsl:if>
        <xsl:if test="$colsConditional = true()"><xsl:text>[ &#160;</xsl:text></xsl:if>
      </td>
      <xsl:call-template name="insertColumns"/>
      <td align="left">
	<xsl:if test="$colsConditional = true()"><xsl:text>&#160; ]</xsl:text></xsl:if>
	<xsl:if test="$colsOptional    = true()"><xsl:text>&#160; }</xsl:text></xsl:if>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="rowgroup" mode="colsMode">
    <xsl:param name="colsOptional"/>
    <xsl:param name="colsConditional"/>
    <xsl:apply-templates select="row | optional | conditional" mode="colsMode">
      <xsl:with-param name="colsOptional" select="$colsOptional"/>
      <xsl:with-param name="colsConditional" select="$colsConditional"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="optional" mode="colsMode">
    <xsl:param name="colsOptional"/>
    <xsl:param name="colsConditional"/>
    <xsl:apply-templates select="row | rowgroup | conditional" mode="colsMode">
      <xsl:with-param name="colsOptional" select="true()"/>
      <xsl:with-param name="colsConditional" select="$colsConditional"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="conditional" mode="colsMode">
    <xsl:param name="colsOptional"/>
    <xsl:param name="colsConditional"/>
    <xsl:apply-templates select="row | rowgroup | optional" mode="colsMode">
      <xsl:with-param name="colsOptional" select="$colsOptional"/>
      <xsl:with-param name="colsConditional" select="true()"/>
    </xsl:apply-templates>
  </xsl:template>
  
  <xsl:template name="insertColumns">
    <xsl:call-template name="insertCol">
      <xsl:with-param name="colID" select="ancestor::cols/@start"/>
    </xsl:call-template>
    <xsl:call-template name="insertCol">
      <xsl:with-param name="colID" select="number(ancestor::cols/@start+1)"/>
    </xsl:call-template>
    <xsl:choose>
      <xsl:when test="number(ancestor::cols/@end) != ancestor::cols/@end">
	<td><xsl:text>&#160;. . .</xsl:text></td>
	<xsl:call-template name="insertCol">
	  <xsl:with-param name="colID" select="ancestor::cols/@end"/>
	</xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
	<xsl:if test="number(ancestor::cols/@end) > number(ancestor::cols/@start+2)">
	  <td><xsl:text>&#160;. . .</xsl:text></td>
	</xsl:if>
	<xsl:if test="number(ancestor::cols/@end) > number(ancestor::cols/@start+1)">
	  <xsl:call-template name="insertCol">
	    <xsl:with-param name="colID" select="ancestor::cols/@end"/>
	  </xsl:call-template>
	</xsl:if>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

  <xsl:template name="insertCol">
    <xsl:param name="colID"/>
    <xsl:message>
      node=<xsl:value-of select="name(.)"/>
    </xsl:message>
    <td>  
      <xsl:text>&#160;</xsl:text>
      <i>
	<xsl:choose>	
	  <xsl:when test="info != '' or status != '' or see != '' or ../../rowgroup/info != ''">
	    <a href="#{generate-id(.)}"><xsl:value-of select="@name"/>(<xsl:value-of select="$colID"/>)</a>
	  </xsl:when>
	  <xsl:otherwise>
	    <xsl:value-of select="@name"/>(<xsl:value-of select="$colID"/>)
	  </xsl:otherwise>
	</xsl:choose>
      </i>
      <xsl:text>&#160;</xsl:text>
    </td>    
  </xsl:template>
  

<!--    *** LINECARD *** -->

  <xsl:template match="linecard">
    <a name="{generate-id(.)}">
      <table border="0" width="100%" style="margin-bottom: 20; ">
	<tr>
	  <th bgcolor="#ddcba6">
	    <h3 style="margin: 10 10 10 15; text-align: left;"> 
	      Line of input 
	    </h3>
	  </th>
	</tr>
	
	<tr><td style="text-align: left; background: #ffebc6; padding: 5 5 5 30; ">	    
	  <table style="border-color: #505087; border-style: solid; border-width: 0; margin-bottom: 10; table-layout: auto; width: 700; ">
	    <tbody>	      	      
	      <tr><td>
		<h3>Syntax:</h3>
		<blockquote>
		  <xsl:apply-templates select="keyword | var | vargroup | list | optional | conditional" mode="syntax"/>
		</blockquote>
	      </td></tr>
	      <tr><td>
		<h3>Description of items:</h3>
		<blockquote>
		  <xsl:apply-templates select="descendant::vargroup | descendant::var | descendant::list"/>
		</blockquote>
	      </td></tr>
	    </tbody>
	  </table>
	</td></tr>
      </table>
    </a>      
  </xsl:template>

  <!--    *** LABEL ***  -->

  <xsl:template match="label">
    <b><xsl:value-of select="."/></b><p/>
  </xsl:template>


  <!--    *** MESSAGE ***  -->

  <xsl:template match="message">
    <pre>
      <xsl:value-of select="."/>
    </pre>
  </xsl:template>


  <!--    *** IF ***  -->

  <xsl:template match="if">
    <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 5">
      <tr><td>
	<b>IF </b>  <xsl:value-of select="@test"/> : <!-- <xsl:apply-templates select="label"/> -->
	<blockquote>
	  <xsl:apply-templates/>
	</blockquote>
      </td></tr>
    </table>
  </xsl:template>


 
 <!--    *** CHOOSE ... ***  -->

  <xsl:template match="choose">
    <table style="border-color: #bb9977; border-style: solid; border-width: 1; margin-bottom: 10; table-layout: auto; width: 100%; padding: 5 5 0 5">
      <tr><td>
	<xsl:apply-templates select="when"/>	
	<xsl:apply-templates select="elsewhen"/>	
	<xsl:apply-templates select="otherwise"/>
      </td></tr>
      <xsl:apply-templates select="message | label" mode="choose"/>
    </table>
  </xsl:template>

  <xsl:template match="when">  
    <b>IF </b> <em><xsl:value-of select="@test"/></em> : <!-- <xsl:apply-templates select="label"/> -->
    <blockquote>
      <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
	<tr><td>
	  <xsl:apply-templates/>
	</td></tr>
      </table>
    </blockquote>
  </xsl:template>

  <xsl:template match="elsewhen">  
    <b>ELSEIF </b> <em><xsl:value-of select="@test"/></em> : <!-- <xsl:apply-templates select="label"/> -->
    <blockquote>
      <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
	<tr><td>
	  <xsl:apply-templates/>	  
	</td></tr>
      </table>
    </blockquote>
  </xsl:template>

  <xsl:template match="otherwise">  
    <b>ELSE </b>
    <blockquote>
      <table style="border-color: #bb9977; border-style: solid; border-width: 3; margin-bottom: 10; table-layout: auto; background-color: #FFddbb; width: 100%; padding: 5 5 0 30">
	<tr><td>	
	  <xsl:apply-templates/>
	</td></tr>
      </table>
    </blockquote>
  </xsl:template>
  
  <!-- *** VARGROUP | DIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup" mode="card_description">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <!--    *** VAR | DIMENSION | LIST ***  -->

  <xsl:template match="var | list | dimension" mode="card_description">
    <!--<xsl:if test="child::node() != ''">-->
    <xsl:if test="info != '' or status != '' or see != ''">
      <xsl:apply-templates select="."/>
    </xsl:if>
  </xsl:template>

  <xsl:template match="var | list | dimension">
    <xsl:if test="name(..) != 'vargroup' and name(..) != 'dimensiongroup'">
      <a name="{generate-id(.)}">
	<a name="{@name}">
	  <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
	    <tr>
	      <xsl:choose>
		<xsl:when test="name(.)='var'">
		  <th align="left" valign="top" width="20%" style="background: #ffff99; padding: 2 2 2 10; ">
		    <xsl:value-of select="@name"/>
		  </th>
		</xsl:when>
		<xsl:when test="name(.)='dimension'">
		  <th width="20%" style="white-space: nowrap; text-align: left; vertical-align: top; background: #ffff99; padding: 2 2 2 10; ">
		    <xsl:value-of select="@name"/>(i), i=<xsl:value-of select="@start"/>,<xsl:value-of select="@end"/>
		  </th>
		</xsl:when>
		<xsl:otherwise>
		  <th width="20%" style="white-space: nowrap; text-align: left; vertical-align: top; background: #ffff99; padding: 2 2 2 10; ">
		    <xsl:value-of select="format"/>
		  </th>
		</xsl:otherwise>
	      </xsl:choose>
	      <td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
		<xsl:value-of select="@type"/>
	      </td>
	    </tr>
	    <xsl:apply-templates select="default"/> 
	    <xsl:apply-templates select="status"/>
	    <xsl:apply-templates select="see"/>
	    <xsl:apply-templates select="info"/>
	  </table>
	  <xsl:call-template name="back_to_top"/>
	</a>
      </a>
    </xsl:if>
  </xsl:template>

  <!-- *** VARGROUP | DIMENSIONGROUP *** -->

  <xsl:template match="vargroup | dimensiongroup">
    <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
      <tr>
	<th align="left" valign="top" width="20%" style="white-space: nowrap; background: #ffff99; padding: 2 2 2 10; ">
	  <xsl:if test="name(.)='vargroup'">
	    <xsl:for-each select="var">
	      <a name="{generate-id(.)}">
		<a name="{@name}">
		  <xsl:value-of select="@name"/><xsl:if test="not(position()=last())">, </xsl:if>
		</a>
	      </a>	
	    </xsl:for-each>	
	  </xsl:if>
	  <xsl:if test="name(.)='dimensiongroup'">
	    <xsl:for-each select="dimension">
	      <a name="{generate-id(.)}">
		<a name="{@name}">
		  <xsl:value-of select="@name"/>(i), 
		  <xsl:if test="position()=last()"> 
		    i=<xsl:value-of select="../@start"/>,<xsl:value-of select="../@end"/>
		  </xsl:if>
		</a>
	      </a>
	    </xsl:for-each>
	  </xsl:if>
	</th>
	<td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	  <xsl:value-of select="@type"/>
	</td>
      </tr>
      <xsl:apply-templates select="default"/> 
      <xsl:apply-templates select="status"/>
      <xsl:apply-templates select="see"/>
      <xsl:apply-templates select="info"/>
    </table>
    <xsl:call-template name="back_to_top"/>
  </xsl:template>

  <!--    *** VAR's elements ***  -->

  <xsl:template match="default">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Default:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:value-of select="."/>
      </td>
    </tr>
  </xsl:template>

  <xsl:template match="status">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>Status:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<xsl:value-of select="."/>
      </td>
    </tr>
  </xsl:template>

 <xsl:template match="see">
    <tr>
      <td style="text-align: right; vertical-align: top; background: #ffffc3; padding: 2 10 2 10; "> <i>See:</i> </td>
      <td style="text-align: left;  vertical-align: top; background: #fff3d9; padding: 2 2 2 5; ">
	<a href="#{normalize-space(.)}"><xsl:value-of select="."/></a>
	<!--<xsl:value-of select="."/>-->
      </td>
    </tr>
 </xsl:template>
 
 <xsl:template match="info">
   <tr><td align="left" valign="top" colspan="2">
     <blockquote>
       <pre><xsl:value-of select="."/></pre>
     </blockquote>
   </td></tr>
 </xsl:template>
 
 
 <!--    *** TABLE ***  -->

 <xsl:template match="table" mode="card_description">
   <xsl:apply-templates select="rows | cols" mode="table"/>   
 </xsl:template>

 <xsl:template match="rows | cols" mode="table">
   <xsl:apply-templates select="col | colgroup | row | rowgroup | optional | conditional" mode="table"/>
 </xsl:template>

 <xsl:template match="optional | conditional" mode="table">
   <xsl:apply-templates select="col | colgroup | row | rowgroup | optional | conditional" mode="table"/>
 </xsl:template>
 
 <xsl:template match="colgroup | rowgroup" mode="table">
   <xsl:if test="info != '' or status != '' or see != ''">
     <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
       <tr>
	 <th width="20%" align="left" valign="top" style="background: #ffff99; padding: 2 2 2 10; ">
	   <xsl:for-each select=".//col | .//row">
	     <a name="{@name}"><a name="{generate-id(.)}">
	       <xsl:value-of select="@name"/>
	     </a></a>
	     <xsl:if test="not(position()=last())">
	       <xsl:text>, </xsl:text>
	     </xsl:if>
	   </xsl:for-each>
	 </th>
	 <td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	   <xsl:value-of select="@type"/>
	 </td>
       </tr>
       <xsl:apply-templates select="default"/> 
       <xsl:apply-templates select="status"/>
       <xsl:apply-templates select="see"/>
       <xsl:apply-templates select="info"/>
     </table>
     <xsl:call-template name="back_to_top"/>
   </xsl:if>
 </xsl:template>
 
 <xsl:template match="col | row" mode="table">
   <!--<xsl:if test="child::node() != ''">-->
   <xsl:if test="info != '' or status != '' or see != ''">
     <table width="100%" style="border-color:   #b5b500; border-style: solid; border-width: 2; margin-bottom: 10; table-layout: auto; background-color: #FFFFFF;">
       <tr>
	 <th align="left" valign="top" width="20%" style="background: #ffff99; padding: 2 2 2 10; ">
	   <a name="{@name}"><a name="{generate-id(.)}">
	     <xsl:value-of select="@name"/>
	   </a></a>
	 </th>
	 <td style="text-align: left; vertical-align: top; background: #ffffc3; padding: 2 2 2 5; ">
	   <xsl:value-of select="@type"/>
	 </td>
       </tr>
       <xsl:apply-templates select="default"/> 
       <xsl:apply-templates select="status"/>
       <xsl:apply-templates select="see"/>
       <xsl:apply-templates select="info"/>
     </table>
     <xsl:call-template name="back_to_top"/>
   </xsl:if>
 </xsl:template>
 
 <!-- *** SECTION *** -->
 <xsl:template match="section">
    <blockquote>
      <a name="{generate-id(.)}">
	<h3><xsl:value-of select="@title"/></h3>
      </a>
      <xsl:apply-templates/>
    </blockquote>
 </xsl:template>

 <xsl:template match="subsection">
   <blockquote>
     <a name="{generate-id(.)}">
       <h4><xsl:value-of select="@title"/></h4>
     </a>
     <xsl:apply-templates/>
   </blockquote>
 </xsl:template>

 <xsl:template match="subsubsection">
   <blockquote>
     <a name="{generate-id(.)}">
       <h5><xsl:value-of select="@title"/></h5>
     </a>
     <xsl:apply-templates/>     
   </blockquote>
 </xsl:template>

 <xsl:template match="paragraph">
   <blockquote>
     <a name="{generate-id(.)}">
       <h6><xsl:value-of select="@title"/></h6>
       <xsl:apply-templates/>
     </a>
   </blockquote>
 </xsl:template>
 
 <xsl:template match="text">
   <blockquote>
     <pre><xsl:value-of select="."/></pre>
   </blockquote>
 </xsl:template>
 
</xsl:stylesheet>
