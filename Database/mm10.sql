-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: mm10
-- ------------------------------------------------------
-- Server version	10.0.35-MariaDB

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `1000genome`
--

DROP TABLE IF EXISTS `1000genome`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `1000genome` (
  `id1000genome` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL,
  `altallelecount` int(11) unsigned DEFAULT NULL,
  `amr_af` float DEFAULT NULL,
  `asn_af` float DEFAULT NULL,
  `afr_af` float DEFAULT NULL,
  `eur_af` float DEFAULT NULL,
  `snpsource` set('LOWCOV','EXOME') DEFAULT NULL,
  PRIMARY KEY (`id1000genome`),
  UNIQUE KEY `unique1000genome` (`chrom`,`start`,`refallele`,`allele`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cadd`
--

DROP TABLE IF EXISTS `cadd`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `cadd` (
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `ref` char(1) NOT NULL,
  `alt` char(1) NOT NULL,
  `rawscore` float NOT NULL,
  `phred` float NOT NULL,
  KEY `cadd_key` (`chrom`,`start`,`ref`,`alt`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `clinvar`
--

DROP TABLE IF EXISTS `clinvar`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `clinvar` (
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `ref` varchar(255) NOT NULL,
  `alt` varchar(255) NOT NULL,
  `rcv` varchar(15) NOT NULL,
  `path` set('drug response','risk factor','not provided','Benign','protective','Likely pathogenic','confers sensitivity','Pathogenic','Uncertain significance','other','Likely benign','association','ClinicalSignificance','Histocompatibility') NOT NULL,
  KEY `chromstart` (`chrom`,`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dgvbp`
--

DROP TABLE IF EXISTS `dgvbp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dgvbp` (
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `depth` int(11) unsigned NOT NULL,
  KEY `dgvbp_chrom` (`chrom`,`start`),
  KEY `dgvbp_cov` (`depth`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `evs`
--

DROP TABLE IF EXISTS `evs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `evs` (
  `idexac` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL,
  `filter` set('PASS','VQSR') DEFAULT NULL,
  `ea_homref` int(11) unsigned DEFAULT NULL,
  `ea_het` int(11) unsigned DEFAULT NULL,
  `ea_homalt` int(11) unsigned DEFAULT NULL,
  `aa_homref` int(11) unsigned DEFAULT NULL,
  `aa_het` int(11) unsigned DEFAULT NULL,
  `aa_homalt` int(11) unsigned DEFAULT NULL,
  PRIMARY KEY (`idexac`),
  UNIQUE KEY `uniqueexac` (`chrom`,`start`,`refallele`,`allele`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `exacGeneScores`
--

DROP TABLE IF EXISTS `exacGeneScores`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exacGeneScores` (
  `transcript` varchar(31) NOT NULL,
  `genesymbol` varchar(150) NOT NULL,
  `chrom` varchar(31) NOT NULL,
  `n_exons` int(10) unsigned NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `end` int(10) unsigned NOT NULL,
  `bp` int(10) unsigned NOT NULL,
  `mu_syn` double NOT NULL,
  `mu_mis` double NOT NULL,
  `mu_lof` double NOT NULL,
  `n_syn` int(10) unsigned NOT NULL,
  `n_mis` int(10) unsigned NOT NULL,
  `n_lof` int(10) unsigned NOT NULL,
  `exp_syn` double NOT NULL,
  `exp_mis` double NOT NULL,
  `exp_lof` double NOT NULL,
  `syn_z` double NOT NULL,
  `mis_z` double NOT NULL,
  `lof_z` double NOT NULL,
  `pLI` double NOT NULL,
  `pRec` double NOT NULL,
  `pNull` double NOT NULL,
  KEY `higenesymbol` (`genesymbol`),
  KEY `chrom` (`chrom`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gap`
--

DROP TABLE IF EXISTS `gap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gap` (
  `bin` smallint(6) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `ix` int(11) NOT NULL,
  `n` char(1) NOT NULL,
  `size` int(10) unsigned NOT NULL,
  `type` varchar(255) NOT NULL,
  `bridge` varchar(255) NOT NULL,
  UNIQUE KEY `chrom_2` (`chrom`(20),`chromStart`),
  KEY `chrom` (`chrom`(20),`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genomicSuperDups`
--

DROP TABLE IF EXISTS `genomicSuperDups`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomicSuperDups` (
  `bin` smallint(6) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `score` int(10) unsigned NOT NULL,
  `strand` char(1) NOT NULL,
  `otherChrom` varchar(255) NOT NULL,
  `otherStart` int(10) unsigned NOT NULL,
  `otherEnd` int(10) unsigned NOT NULL,
  `otherSize` int(10) unsigned NOT NULL,
  `uid` int(10) unsigned NOT NULL,
  `posBasesHit` int(10) unsigned NOT NULL,
  `testResult` varchar(255) NOT NULL,
  `verdict` varchar(255) NOT NULL,
  `chits` varchar(255) NOT NULL,
  `ccov` varchar(255) NOT NULL,
  `alignfile` varchar(255) NOT NULL,
  `alignL` int(10) unsigned NOT NULL,
  `indelN` int(10) unsigned NOT NULL,
  `indelS` int(10) unsigned NOT NULL,
  `alignB` int(10) unsigned NOT NULL,
  `matchB` int(10) unsigned NOT NULL,
  `mismatchB` int(10) unsigned NOT NULL,
  `transitionsB` int(10) unsigned NOT NULL,
  `transversionsB` int(10) unsigned NOT NULL,
  `fracMatch` float NOT NULL,
  `fracMatchIndel` float NOT NULL,
  `jcK` float NOT NULL,
  `k2K` float NOT NULL,
  KEY `name` (`name`(32)),
  KEY `chrom` (`chrom`(8),`bin`),
  KEY `chrom_2` (`chrom`(8),`chromStart`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `kaviar`
--

DROP TABLE IF EXISTS `kaviar`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `kaviar` (
  `idkaviar` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL,
  `af` float DEFAULT NULL,
  `ac` int(11) unsigned DEFAULT NULL,
  `an` int(11) unsigned DEFAULT NULL,
  `ds` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`idkaviar`),
  UNIQUE KEY `uniquekaviar` (`chrom`,`start`,`refallele`,`allele`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `kgXref`
--

DROP TABLE IF EXISTS `kgXref`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `kgXref` (
  `kgID` varchar(255) NOT NULL,
  `mRNA` varchar(255) NOT NULL,
  `spID` varchar(255) NOT NULL,
  `spDisplayID` varchar(255) NOT NULL,
  `geneSymbol` varchar(255) NOT NULL,
  `refseq` varchar(255) NOT NULL,
  `protAcc` varchar(255) NOT NULL,
  `description` longblob NOT NULL,
  `rfamAcc` varchar(255) NOT NULL,
  `tRnaName` varchar(255) NOT NULL,
  KEY `kgID` (`kgID`),
  KEY `mRNA` (`mRNA`),
  KEY `spID` (`spID`),
  KEY `spDisplayID` (`spDisplayID`),
  KEY `geneSymbol` (`geneSymbol`),
  KEY `refseq` (`refseq`),
  KEY `protAcc` (`protAcc`),
  KEY `rfamAcc` (`rfamAcc`),
  KEY `tRnaName` (`tRnaName`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownGene`
--

DROP TABLE IF EXISTS `knownGene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGene` (
  `name` varchar(255) NOT NULL DEFAULT '',
  `chrom` varchar(255) NOT NULL DEFAULT '',
  `strand` char(1) NOT NULL DEFAULT '',
  `txStart` int(10) unsigned NOT NULL DEFAULT '0',
  `txEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `cdsStart` int(10) unsigned NOT NULL DEFAULT '0',
  `cdsEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `exonCount` int(10) unsigned NOT NULL DEFAULT '0',
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `proteinID` varchar(40) NOT NULL DEFAULT '',
  `alignID` varchar(255) NOT NULL DEFAULT '',
  KEY `name` (`name`),
  KEY `chrom` (`chrom`(16),`txStart`),
  KEY `chrom_2` (`chrom`(16),`txEnd`),
  KEY `protein` (`proteinID`(16)),
  KEY `align` (`alignID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownGenePep`
--

DROP TABLE IF EXISTS `knownGenePep`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGenePep` (
  `name` varchar(255) NOT NULL,
  `seq` longblob NOT NULL,
  PRIMARY KEY (`name`(64))
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary table structure for view `knownGeneSymbol`
--

DROP TABLE IF EXISTS `knownGeneSymbol`;
/*!50001 DROP VIEW IF EXISTS `knownGeneSymbol`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `knownGeneSymbol` (
  `name` tinyint NOT NULL,
  `chrom` tinyint NOT NULL,
  `strand` tinyint NOT NULL,
  `txStart` tinyint NOT NULL,
  `txEnd` tinyint NOT NULL,
  `cdsStart` tinyint NOT NULL,
  `cdsEnd` tinyint NOT NULL,
  `exonCount` tinyint NOT NULL,
  `exonStarts` tinyint NOT NULL,
  `exonEnds` tinyint NOT NULL,
  `proteinID` tinyint NOT NULL,
  `alignID` tinyint NOT NULL,
  `geneSymbol` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `knownGene_cds`
--

DROP TABLE IF EXISTS `knownGene_cds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGene_cds` (
  `name` varchar(255) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pph3`
--

DROP TABLE IF EXISTS `pph3`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pph3` (
  `chrom` varchar(45) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `ref` char(1) DEFAULT NULL,
  `alt` char(1) DEFAULT NULL,
  `transcript` varchar(45) DEFAULT NULL,
  `str` char(1) DEFAULT NULL,
  `gene` varchar(45) DEFAULT NULL,
  `refs_acc` varchar(45) DEFAULT NULL,
  `cdnpos` tinyint(4) DEFAULT NULL,
  `frame` tinyint(4) DEFAULT NULL,
  `nt1` char(1) DEFAULT NULL,
  `nt2` char(1) DEFAULT NULL,
  `rsid` varchar(45) DEFAULT NULL,
  `acc` varchar(45) DEFAULT NULL,
  `pos` int(11) DEFAULT NULL,
  `aa1` char(1) DEFAULT NULL,
  `aa2` char(1) DEFAULT NULL,
  `hdiv_prediction` varchar(45) DEFAULT NULL,
  `hdiv_class` varchar(45) DEFAULT NULL,
  `hdiv_prob` float DEFAULT NULL,
  `hdiv_FPR` float DEFAULT NULL,
  `hdiv_TPR` float DEFAULT NULL,
  `hdiv_FDR` float DEFAULT NULL,
  `hvar_prediction` varchar(45) DEFAULT NULL,
  `hvar_class` varchar(45) DEFAULT NULL,
  `hvar_prob` float DEFAULT NULL,
  `hvar_FPR` float DEFAULT NULL,
  `hvar_TPR` float DEFAULT NULL,
  `hvar_FDR` float DEFAULT NULL,
  KEY `pph3_key` (`chrom`,`start`,`ref`,`alt`),
  KEY `pph3_transcsript` (`transcript`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sift`
--

DROP TABLE IF EXISTS `sift`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sift` (
  `chrom` varchar(45) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `snp` varchar(45) DEFAULT NULL,
  `ref` char(1) DEFAULT NULL,
  `alt` char(1) DEFAULT NULL,
  `score` float DEFAULT NULL,
  `median` float DEFAULT NULL,
  KEY `sift_key` (`chrom`,`start`,`ref`,`alt`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snp137`
--

DROP TABLE IF EXISTS `snp137`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp137` (
  `bin` smallint(5) unsigned NOT NULL,
  `chrom` varchar(31) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `name` varchar(15) NOT NULL,
  `score` smallint(5) unsigned NOT NULL,
  `strand` enum('+','-') NOT NULL,
  `refNCBI` blob NOT NULL,
  `refUCSC` blob NOT NULL,
  `observed` varchar(255) NOT NULL,
  `molType` enum('genomic','cDNA','mito') NOT NULL,
  `class` enum('single','in-del','named','mixed','mnp','insertion','deletion') NOT NULL,
  `valid` set('unknown','by-cluster','by-frequency') NOT NULL,
  `avHet` float NOT NULL,
  `avHetSE` float NOT NULL,
  `func` set('unknown','coding-synon','intron','near-gene-3','near-gene-5','ncRNA','nonsense','missense','stop-loss','frameshift','cds-indel','untranslated-3','untranslated-5','splice-3','splice-5') NOT NULL,
  `locType` enum('range','exact','between','rangeDeletion') NOT NULL,
  `weight` int(10) unsigned NOT NULL,
  `exceptions` set('RefAlleleMismatch','DuplicateObserved','FlankMismatchGenomeShorter','SingleClassLongerSpan','SingleClassZeroSpan','SingleClassTriAllelic','SingleClassQuadAllelic','ObservedWrongFormat','ObservedTooLong','ObservedContainsIupac','ObservedMismatch','MultipleAlignments','SingleAlleleFreq','InconsistentAlleles') NOT NULL,
  `submitterCount` smallint(5) unsigned NOT NULL,
  `submitters` longblob NOT NULL,
  `alleleFreqCount` smallint(5) unsigned NOT NULL,
  `alleles` longblob NOT NULL,
  `alleleNs` longblob NOT NULL,
  `alleleFreqs` longblob NOT NULL,
  `bitfields` set('maf-5-some-pop','maf-5-all-pops','genotype-conflict','rs-cluster-nonoverlapping-alleles','observed-mismatch') NOT NULL,
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snp142`
--

DROP TABLE IF EXISTS `snp142`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp142` (
  `bin` smallint(5) unsigned NOT NULL,
  `chrom` varchar(31) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `name` varchar(15) NOT NULL,
  `score` smallint(5) unsigned NOT NULL,
  `strand` enum('+','-') NOT NULL,
  `refNCBI` blob NOT NULL,
  `refUCSC` blob NOT NULL,
  `observed` varchar(255) NOT NULL,
  `molType` enum('genomic','cDNA','mito') NOT NULL,
  `class` enum('single','in-del','mnp','insertion','deletion') NOT NULL,
  `valid` set('unknown','by-cluster','by-frequency') NOT NULL,
  `avHet` float NOT NULL,
  `avHetSE` float NOT NULL,
  `func` set('unknown','coding-synon','intron','near-gene-3','near-gene-5','ncRNA','nonsense','missense','stop-loss','frameshift','cds-indel','untranslated-3','untranslated-5','splice-3','splice-5') NOT NULL,
  `locType` enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion','fuzzy') NOT NULL,
  `weight` int(10) unsigned NOT NULL,
  `exceptions` set('RefAlleleMismatch','DuplicateObserved','MixedObserved','FlankMismatchGenomeLonger','FlankMismatchGenomeEqual','FlankMismatchGenomeShorter','SingleClassLongerSpan','SingleClassZeroSpan','SingleClassTriAllelic','SingleClassQuadAllelic','ObservedWrongFormat','ObservedTooLong','ObservedContainsIupac','ObservedMismatch','MultipleAlignments','SingleAlleleFreq','InconsistentAlleles') NOT NULL,
  `submitterCount` smallint(5) unsigned NOT NULL,
  `submitters` longblob NOT NULL,
  `alleleFreqCount` smallint(5) unsigned NOT NULL,
  `alleles` longblob NOT NULL,
  `alleleNs` longblob NOT NULL,
  `alleleFreqs` longblob NOT NULL,
  `bitfields` set('maf-5-some-pop','maf-5-all-pops','genotype-conflict','rs-cluster-nonoverlapping-alleles','observed-mismatch') NOT NULL,
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Final view structure for view `knownGeneSymbol`
--

/*!50001 DROP TABLE IF EXISTS `knownGeneSymbol`*/;
/*!50001 DROP VIEW IF EXISTS `knownGeneSymbol`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`wieland`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `knownGeneSymbol` AS select `kg`.`name` AS `name`,`kg`.`chrom` AS `chrom`,`kg`.`strand` AS `strand`,`kg`.`txStart` AS `txStart`,`kg`.`txEnd` AS `txEnd`,`kg`.`cdsStart` AS `cdsStart`,`kg`.`cdsEnd` AS `cdsEnd`,`kg`.`exonCount` AS `exonCount`,`kg`.`exonStarts` AS `exonStarts`,`kg`.`exonEnds` AS `exonEnds`,`kg`.`proteinID` AS `proteinID`,`kg`.`alignID` AS `alignID`,`x`.`geneSymbol` AS `geneSymbol` from (`knownGene` `kg` join `kgXref` `x` on((`x`.`kgID` = `kg`.`name`))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2019-11-25 13:57:14
