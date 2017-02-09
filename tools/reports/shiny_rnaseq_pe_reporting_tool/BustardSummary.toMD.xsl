<?xml version="1.0"?> 
<!--
This file is an XSL program to transform Bustard XML files into Latex.
It is thought to be run in conjuntion with the xmlparser.pl to generate
the Markdown tables included in the reports.
Who:  Sergi Sayols
When: 27-01-2014
-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"> 
<xsl:output method="text"/>
<xsl:template match="/">

<!-- Chip summary report -->
<xsl:for-each select="BustardSummary/ChipSummary"> 
Desc | Value
:--: | :---:
Machine | <xsl:value-of select="Machine"/>
Run Folder | <xsl:variable name="RF" select="RunFolder"/> <xsl:value-of select="translate($RF,'_','-')" />
Chip ID | <xsl:value-of select="concat(ChipID,' &#10;')"/>
</xsl:for-each> 

<!-- Samples summary report -->
<xsl:if test="string(BustardSummary/Samples)">
Lane | Barcode | Sample | Species
:--: | :-----: | :----: | :-----:
<xsl:for-each select="BustardSummary/Samples/Lane">
<xsl:value-of select="laneNumber"/> | <xsl:value-of select="barcode"/> | <xsl:value-of select="sampleId"/> | <xsl:value-of select="concat(species,' &#10;')"/>
</xsl:for-each>
</xsl:if>

<!-- Lane results summary report -->
<xsl:for-each select="BustardSummary/LaneResultsSummary/Read">
<xsl:variable name="numReads" select="count(../Read)"/>
<xsl:variable name="r" select="readNumber" />
<xsl:variable name="Caption" select="concat('Lane Results Summary: Read ',$r)" />
Lane | Lane Yield (kbases) | Clusters (raw) | Clusters (PF) | First Cycle Int (PF) | \% intensity after 20 cycles (PF) | \% PF Clusters
:--: | :-----------------: | :------------: | :-----------: | :------------------: | :--------------------------: | :------------:
<xsl:variable name="clusterCountRawMean" select="sum(Lane/clusterCountRaw/mean)"/>
<xsl:variable name="clusterCountPFMean"  select="sum(Lane/clusterCountPF/mean)"/>
<xsl:variable name="oneSigMean"			 select="sum(Lane/oneSig/mean)"/>
<xsl:variable name="signal20AsPctOf1Mean"  select="sum(Lane/signal20AsPctOf1/mean)"/>
<xsl:variable name="percentClustersPFMean" select="sum(Lane/percentClustersPF/mean)"/>
<xsl:variable name="numLanes"		 	 select="count(Lane/laneYield)"/>
<xsl:for-each select="Lane">
<xsl:if test="string(laneYield)">
<xsl:value-of select="laneNumber"/>  | <xsl:value-of select="laneYield"/> | <xsl:value-of select="clusterCountRaw/mean"/> $\pm$ <xsl:value-of select="clusterCountRaw/stdev"/> | <xsl:value-of select="clusterCountPF/mean"/> $\pm$ <xsl:value-of select="clusterCountPF/stdev"/> | <xsl:value-of select="oneSig/mean"/> $\pm$ <xsl:value-of select="oneSig/stdev"/> | <xsl:value-of select="signal20AsPctOf1/mean"/> $\pm$ <xsl:value-of select="signal20AsPctOf1/stdev"/> | <xsl:value-of select="percentClustersPF/mean"/> $\pm$ <xsl:value-of select="concat(percentClustersPF/stdev,' &#10;')"/>
</xsl:if>
</xsl:for-each>
_ | **Average** | <xsl:value-of select="round($clusterCountRawMean div $numLanes)"/> | <xsl:value-of select="round($clusterCountPFMean div $numLanes)"/> | <xsl:value-of select="round($oneSigMean div $numLanes)"/> | <xsl:value-of select="round($signal20AsPctOf1Mean div $numLanes * 100) div 100"/> | <xsl:value-of select="concat(round($percentClustersPFMean div $numLanes * 100) div 100,' &#10;')"/>
</xsl:for-each>

<!-- Expanded lane summary report -->
<xsl:for-each select="BustardSummary/ExpandedLaneSummary/Read">
	<xsl:variable name="numReads" select="count(../Read)"/>
	<xsl:variable name="r" select="readNumber" />
	<xsl:variable name="Caption" select="concat('Expanded Lane Summary: Read ',$r)" />
Lane | Clusters(tile $\mu$,raw) | \% Phas | \% Prephas | \% Retained(raw) | Cyc2-4 $\mu$ Int(raw,PF) | Cyc2-10 $\mu$ \% Loss(filt,PF) | Cyc10-20 $\mu$ \% Loss(filt,PF)
:--: | :----------------------: | :-----: | :--------: | :--------------: | :------------------------: | :----------------------------: | :-----------------------------:
<xsl:for-each select="Lane">
<xsl:if test="string(phasingApplied)">
<xsl:value-of select="laneNumber"/> | <xsl:value-of select="clusterCountRaw/mean"/> | <xsl:value-of select="phasingApplied"/> | <xsl:value-of select="prephasingApplied"/> | <xsl:value-of select="percentClustersPF/mean"/> | <xsl:value-of select="signalAverage2to4/mean"/> $\pm$ <xsl:value-of select="signalAverage2to4/stdev"/> | <xsl:value-of select="signalLoss2to10/mean"/> $\pm$ <xsl:value-of select="signalLoss2to10/stdev"/> | <xsl:value-of select="signalLoss10to20/mean"/> $\pm$ <xsl:value-of select="concat(signalLoss10to20/stdev,' &#10;')"/>
</xsl:if>
</xsl:for-each>
</xsl:for-each>

</xsl:template>
</xsl:stylesheet>
