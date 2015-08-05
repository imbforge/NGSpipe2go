<?xml version="1.0"?> 
<!--
This file is an XSL program to transform Bustard XML files into Latex.
It is thought to be run in conjuntion with the xmlparser.pl to generate
the Latex tables included in the reports.
Who:  Sergi Sayols
When: 27-01-2014
-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0"> 
<xsl:output method="text"/>

<xsl:template match="/">

\subsection{SAV}

<!-- Chip summary report -->
	\begin{table}[h!]
	\tiny
	\begin{tabular}{|c|c|}\hline
	<xsl:for-each select="BustardSummary/ChipSummary"> 
		Machine <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="Machine"/> \\\hline
		Run Folder <xsl:value-of select=" '&amp; '"/>
		<xsl:variable name="RF" select="RunFolder"/> 
		<xsl:value-of select="translate($RF,'_','-')" /> \\\hline
		Chip ID <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="ChipID"/> \\\hline
	</xsl:for-each> 
	\end{tabular}
	\caption{Chip Summary}\label{Chip Summary}
	\end{table}

<!-- Samples summary report -->
<xsl:if test="string(BustardSummary/Samples)">
	\begin{table}[h!]
	\tiny
	\begin{tabular}{|c|c|c|c|}\hline
		\textbf{Lane} <xsl:value-of select=" '&amp; '"/>
		\textbf{Barcode} <xsl:value-of select=" '&amp; '"/>
		\textbf{Sample} <xsl:value-of select=" '&amp; '"/>
		\textbf{Species} <xsl:value-of select=" '&amp; '"/> \\\hline
	<xsl:for-each select="BustardSummary/Samples/Lane">
		 <xsl:value-of select="laneNumber"/> <xsl:value-of select=" '&amp; '"/>
		 <xsl:value-of select="barcode"/> <xsl:value-of select=" '&amp; '"/>
		 <xsl:value-of select="sampleId"/> <xsl:value-of select=" '&amp; '"/>
		 <xsl:value-of select="species"/> \\\hline
	</xsl:for-each>
	\end{tabular}
	\caption{Samples Summary}\label{Samples Summary}
	\end{table}
</xsl:if>

<!-- Lane results summary report -->
<xsl:for-each select="BustardSummary/LaneResultsSummary/Read">
	<xsl:variable name="numReads" select="count(../Read)"/>
	<xsl:variable name="r" select="readNumber" />
	<xsl:variable name="Caption" select="concat('Lane Results Summary: Read ',$r)" />
	\begin{table}[h!]
	\tiny
	\begin{tabular}{|c|c|c|c|c|c|c|}\hline
		\textbf{Lane} <xsl:value-of select=" '&amp; '"/>
		\textbf{Lane Yield (kbases)} <xsl:value-of select=" '&amp; '"/>
		\textbf{Clusters (raw)} <xsl:value-of select=" '&amp; '"/>
		\textbf{Clusters (PF)} <xsl:value-of select=" '&amp; '"/>
		\textbf{First Cycle Int (PF)} <xsl:value-of select=" '&amp; '"/>
		\textbf{\% intensity after 20 cycles (PF)} <xsl:value-of select=" '&amp; '"/>
		\textbf{\% PF Clusters} \\\hline

		<xsl:variable name="clusterCountRawMean" select="sum(Lane/clusterCountRaw/mean)"/>
		<xsl:variable name="clusterCountPFMean"  select="sum(Lane/clusterCountPF/mean)"/>
		<xsl:variable name="oneSigMean"			 select="sum(Lane/oneSig/mean)"/>
		<xsl:variable name="signal20AsPctOf1Mean"  select="sum(Lane/signal20AsPctOf1/mean)"/>
		<xsl:variable name="percentClustersPFMean" select="sum(Lane/percentClustersPF/mean)"/>
		<xsl:variable name="numLanes"		 	 select="count(Lane/laneYield)"/>

		<xsl:for-each select="Lane">
			<xsl:if test="string(laneYield)">
				<xsl:value-of select="laneNumber"/>  <xsl:value-of select=" '&amp; '"/>
				<xsl:value-of select="laneYield"/> <xsl:value-of select=" '&amp; '"/>
				<xsl:value-of select="clusterCountRaw/mean"/> $\pm$ <xsl:value-of select="clusterCountRaw/stdev"/> <xsl:value-of select=" '&amp; '"/>
				<xsl:value-of select="clusterCountPF/mean"/> $\pm$ <xsl:value-of select="clusterCountPF/stdev"/> <xsl:value-of select=" '&amp; '"/>
				<xsl:value-of select="oneSig/mean"/> $\pm$ <xsl:value-of select="oneSig/stdev"/> <xsl:value-of select=" '&amp; '"/>
				<xsl:value-of select="signal20AsPctOf1/mean"/> $\pm$ <xsl:value-of select="signal20AsPctOf1/stdev"/> <xsl:value-of select=" '&amp; '"/>
				<xsl:value-of select="percentClustersPF/mean"/> $\pm$ <xsl:value-of select="percentClustersPF/stdev"/> \\\hline
			</xsl:if>
		</xsl:for-each>

		<xsl:value-of select=" '&amp; '"/> \textbf{Average} <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="round($clusterCountRawMean div $numLanes)"/> <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="round($clusterCountPFMean div $numLanes)"/> <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="round($oneSigMean div $numLanes)"/> <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="round($signal20AsPctOf1Mean div $numLanes * 100) div 100"/> <xsl:value-of select=" '&amp; '"/>
		<xsl:value-of select="round($percentClustersPFMean div $numLanes * 100) div 100"/> \\\hline
	\end{tabular}
	\caption{<xsl:value-of select="$Caption" />}\label{<xsl:value-of select="$Caption" />}
	\end{table}
</xsl:for-each>


<!-- Expanded lane summary report -->
<xsl:for-each select="BustardSummary/ExpandedLaneSummary/Read">
	<xsl:variable name="numReads" select="count(../Read)"/>
	<xsl:variable name="r" select="readNumber" />
	<xsl:variable name="Caption" select="concat('Expanded Lane Summary: Read ',$r)" />
	\begin{table}[h!]
	\tiny
	\begin{tabular}{|c|c|c|c|c|c|c|c|}\hline
	\textbf{Lane} <xsl:value-of select=" '&amp; '"/>
	\textbf{Clusters(tile $\mu$,raw)} <xsl:value-of select=" '&amp; '"/>
	\textbf{\% Phas} <xsl:value-of select=" '&amp; '"/>
	\textbf{\% Prephas} <xsl:value-of select=" '&amp; '"/>
	\textbf{\% Retained(raw)} <xsl:value-of select=" '&amp; '"/>
	\textbf{Cyc2-4 $\mu$ Int(raw,PF)} <xsl:value-of select=" '&amp; '"/>
	\textbf{Cyc2-10 $\mu$ \% Loss(filt,PF)} <xsl:value-of select=" '&amp; '"/>
	\textbf{Cyc10-20 $\mu$ \% Loss(filt,PF)} \\\hline

	<xsl:for-each select="Lane">
		<xsl:if test="string(phasingApplied)">
			<xsl:value-of select="laneNumber"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="clusterCountRaw/mean"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="phasingApplied"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="prephasingApplied"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="percentClustersPF/mean"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="signalAverage2to4/mean"/> $\pm$ <xsl:value-of select="signalAverage2to4/stdev"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="signalLoss2to10/mean"/> $\pm$ <xsl:value-of select="signalLoss2to10/stdev"/> <xsl:value-of select=" '&amp; '"/>
			<xsl:value-of select="signalLoss10to20/mean"/> $\pm$ <xsl:value-of select="signalLoss10to20/stdev"/> \\\hline
		</xsl:if>

	</xsl:for-each>
	\end{tabular}
	\caption{<xsl:value-of select="$Caption" />}\label{<xsl:value-of select="$Caption" />}
	\end{table}
</xsl:for-each>

</xsl:template>

</xsl:stylesheet>
