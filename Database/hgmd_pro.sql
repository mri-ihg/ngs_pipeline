-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: hgmd_pro
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
-- Table structure for table `allgenes`
--

DROP TABLE IF EXISTS `allgenes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `allgenes` (
  `disease` text COLLATE utf8_unicode_ci,
  `gene_id` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `chrom` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `genename` varchar(200) COLLATE utf8_unicode_ci DEFAULT NULL,
  `refseq` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `altsymbol` text COLLATE utf8_unicode_ci,
  `altname` text COLLATE utf8_unicode_ci,
  `gdbid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `omimid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `entrezID` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `hgncID` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `svar` int(11) DEFAULT NULL,
  `mut` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `poly` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `ftv` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `go_terms_acc` text COLLATE utf8_unicode_ci,
  `go_terms_name` text COLLATE utf8_unicode_ci,
  `mut_total` int(11) DEFAULT NULL,
  `new_mut_total` int(11) DEFAULT NULL,
  `gene_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`),
  UNIQUE KEY `gene_id` (`gene_id`),
  UNIQUE KEY `refseq` (`refseq`),
  UNIQUE KEY `gdbid` (`gdbid`),
  UNIQUE KEY `entrezID` (`entrezID`),
  UNIQUE KEY `hgncID` (`hgncID`),
  FULLTEXT KEY `ALLFIELDS` (`disease`,`gene`,`altsymbol`,`altname`,`chrom`,`genename`,`gdbid`,`omimid`,`entrezID`,`hgncID`,`refseq`,`comments`,`go_terms_acc`,`go_terms_name`),
  FULLTEXT KEY `GENEDATA` (`gene`,`altsymbol`,`altname`,`genename`),
  FULLTEXT KEY `OTHERDATA` (`hgncID`,`omimid`,`gdbid`,`entrezID`),
  FULLTEXT KEY `GOTERMS` (`go_terms_acc`,`go_terms_name`),
  FULLTEXT KEY `DISEASE` (`disease`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `allmut`
--

DROP TABLE IF EXISTS `allmut`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `allmut` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `chrom` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `genename` varchar(200) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gdbid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `omimid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `amino` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `deletion` varchar(65) COLLATE utf8_unicode_ci DEFAULT NULL,
  `insertion` varchar(65) COLLATE utf8_unicode_ci DEFAULT NULL,
  `codon` int(11) DEFAULT NULL,
  `codonAff` int(11) DEFAULT NULL,
  `descr` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `hgvs` varchar(60) COLLATE utf8_unicode_ci DEFAULT NULL,
  `hgvsAll` varchar(120) COLLATE utf8_unicode_ci DEFAULT NULL,
  `dbsnp` varchar(12) COLLATE utf8_unicode_ci DEFAULT NULL,
  `chromosome` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `startCoord` int(11) DEFAULT NULL,
  `endCoord` int(11) DEFAULT NULL,
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `dmsupport` int(11) DEFAULT NULL,
  `rankscore` double DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `allname` varchar(200) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` char(4) COLLATE utf8_unicode_ci DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `pmidAll` text COLLATE utf8_unicode_ci,
  `reftag` char(3) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  `base` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`acc_num`),
  KEY `hgvs` (`hgvs`),
  KEY `pmid` (`pmid`),
  KEY `gene` (`gene`),
  FULLTEXT KEY `ALLFIELDS` (`disease`,`gene`,`chrom`,`genename`,`gdbid`,`omimid`,`author`,`fullname`,`allname`,`vol`,`page`,`year`,`pmid`,`pmidAll`,`comments`,`acc_num`),
  FULLTEXT KEY `REFFIELDS` (`gene`,`author`,`fullname`,`allname`,`year`,`pmid`,`pmidAll`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `altname`
--

DROP TABLE IF EXISTS `altname`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `altname` (
  `altname_id` int(11) NOT NULL DEFAULT '0',
  `gene_id` int(11) DEFAULT NULL,
  `altname` varchar(100) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`altname_id`),
  UNIQUE KEY `gene_id` (`gene_id`,`altname`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `altsym`
--

DROP TABLE IF EXISTS `altsym`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `altsym` (
  `alt_id` int(11) NOT NULL DEFAULT '0',
  `gene_id` int(11) DEFAULT NULL,
  `altsymbol` varchar(35) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`alt_id`),
  UNIQUE KEY `gene_id` (`gene_id`,`altsymbol`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `amino`
--

DROP TABLE IF EXISTS `amino`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `amino` (
  `code3` char(4) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `code1` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `full` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `polarity` char(3) COLLATE utf8_unicode_ci DEFAULT NULL,
  `pH` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `weight` int(11) DEFAULT NULL,
  `hydrophob` decimal(2,1) DEFAULT NULL,
  `hydrophil` decimal(2,1) DEFAULT NULL,
  `sec` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`code3`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `aminoseq`
--

DROP TABLE IF EXISTS `aminoseq`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `aminoseq` (
  `gene_id` int(11) NOT NULL DEFAULT '0',
  `gene` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `txid` int(11) DEFAULT NULL,
  `species` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `sp_name` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `prot_seq` text COLLATE utf8_unicode_ci,
  `prot_acc` varchar(25) COLLATE utf8_unicode_ci NOT NULL,
  `prot_ver` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `curated` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`prot_acc`,`gene_id`),
  KEY `gene` (`gene`),
  KEY `txid` (`txid`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `amplet`
--

DROP TABLE IF EXISTS `amplet`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `amplet` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `chrom` varchar(3) COLLATE utf8_unicode_ci DEFAULT NULL,
  `amplet` varchar(50) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `norcopy` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `patcopy` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `location` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`amplet`,`norcopy`,`patcopy`,`location`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `complex`
--

DROP TABLE IF EXISTS `complex`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `complex` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `descr` varchar(75) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`descr`,`pmid`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dbnsfp2`
--

DROP TABLE IF EXISTS `dbnsfp2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dbnsfp2` (
  `id` int(11) NOT NULL DEFAULT '0',
  `acc_num` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `uniprot` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `polyphen` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `LRT` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutationTaster` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutationAssessor` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `FATHMM` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `phyloP` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  `GERP_RS` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dbnsfp3`
--

DROP TABLE IF EXISTS `dbnsfp3`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dbnsfp3` (
  `id` int(11) NOT NULL DEFAULT '0',
  `acc_num` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `SIFT` varchar(120) COLLATE utf8_unicode_ci DEFAULT NULL,
  `LRT` varchar(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutationTaster` varchar(67) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutationAssessor` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `FATHMM` varchar(120) COLLATE utf8_unicode_ci DEFAULT NULL,
  `PROVEAN` varchar(120) COLLATE utf8_unicode_ci DEFAULT NULL,
  `MCAP` varchar(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fathmm_MKL` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `MetaSVM_pred` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `MetaLR_pred` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `ExAC_AC` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `ExAC_AF` varchar(9) COLLATE utf8_unicode_ci DEFAULT NULL,
  `1000Gp3_AC` varchar(4) COLLATE utf8_unicode_ci DEFAULT NULL,
  `1000Gp3_AF` varchar(21) COLLATE utf8_unicode_ci DEFAULT NULL,
  `Interpro` text COLLATE utf8_unicode_ci,
  `phastCons100` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `phyloP100` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `phastCons20` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `phyloP20` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `GERP_RS` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dbsnp`
--

DROP TABLE IF EXISTS `dbsnp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dbsnp` (
  `hgmd_acc` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `dbsnp_id` varchar(12) COLLATE utf8_unicode_ci DEFAULT NULL,
  `het` decimal(2,2) DEFAULT NULL,
  `MAFsource` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `MAFfreq` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  PRIMARY KEY (`hgmd_acc`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `deletion`
--

DROP TABLE IF EXISTS `deletion`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `deletion` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `deletion` varchar(65) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `codon` int(11) NOT NULL DEFAULT '0',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`deletion`,`codon`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `delins`
--

DROP TABLE IF EXISTS `delins`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `delins` (
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `delins` varchar(65) COLLATE utf8_unicode_ci DEFAULT NULL,
  `deleted` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `inserted` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `del_length` int(2) DEFAULT NULL,
  `ins_length` int(2) DEFAULT NULL,
  `base` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `extrarefs`
--

DROP TABLE IF EXISTS `extrarefs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `extrarefs` (
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `disease` varchar(125) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `gene` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `title` text COLLATE utf8_unicode_ci,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `allname` varchar(200) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `issue` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `year` char(4) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `support` enum('1','0','-1') COLLATE utf8_unicode_ci DEFAULT NULL,
  `reftag` char(4) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`acc_num`,`author`,`fullname`,`vol`,`page`,`year`,`reftag`,`disease`),
  FULLTEXT KEY `REFFIELDS` (`gene`,`author`,`fullname`,`allname`,`year`,`pmid`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `faq`
--

DROP TABLE IF EXISTS `faq`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `faq` (
  `id` int(10) unsigned NOT NULL,
  `question` text COLLATE utf8_unicode_ci,
  `answer` text COLLATE utf8_unicode_ci,
  `owner` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `type` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `func_anotat`
--

DROP TABLE IF EXISTS `func_anotat`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `func_anotat` (
  `acc_num` varchar(30) COLLATE utf8_unicode_ci NOT NULL,
  `site_id` int(11) NOT NULL,
  `pmid` int(11) NOT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `direction` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `wildtype` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutant` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`acc_num`,`site_id`,`pmid`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `func_sites`
--

DROP TABLE IF EXISTS `func_sites`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `func_sites` (
  `id` int(11) NOT NULL,
  `abr` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `name` varchar(200) COLLATE utf8_unicode_ci NOT NULL,
  `type` varchar(60) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` varchar(250) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `new_index` (`type`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene2refseq`
--

DROP TABLE IF EXISTS `gene2refseq`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene2refseq` (
  `hgmdID` int(11) NOT NULL DEFAULT '0',
  `refcore` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `refversion` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `offset` int(11) DEFAULT NULL,
  `first` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`hgmdID`),
  UNIQUE KEY `refcore` (`refcore`,`refversion`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genefam`
--

DROP TABLE IF EXISTS `genefam`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genefam` (
  `url` varchar(100) COLLATE utf8_unicode_ci DEFAULT NULL,
  `genetag` varchar(25) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `genefamily` varchar(150) COLLATE utf8_unicode_ci DEFAULT NULL,
  `symbol` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `hgncID` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`genetag`,`hgncID`),
  KEY `symbol` (`symbol`),
  KEY `hgncID` (`hgncID`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `goterms`
--

DROP TABLE IF EXISTS `goterms`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `goterms` (
  `gene_id` int(11) NOT NULL DEFAULT '0',
  `hgmdsym` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `go_term_name` varchar(255) COLLATE utf8_unicode_ci DEFAULT NULL,
  `go_term_type` varchar(55) COLLATE utf8_unicode_ci DEFAULT NULL,
  `go_term_acc` varchar(100) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `go_term_evid` varchar(8) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`gene_id`,`go_term_acc`,`go_term_evid`),
  FULLTEXT KEY `GOFIELDS` (`hgmdsym`,`go_term_name`,`go_term_type`,`go_term_acc`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `grosdel`
--

DROP TABLE IF EXISTS `grosdel`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `grosdel` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `descr` varchar(75) COLLATE utf8_unicode_ci DEFAULT NULL,
  `cdna` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`acc_num`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `grosins`
--

DROP TABLE IF EXISTS `grosins`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `grosins` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `type` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `descr` varchar(75) COLLATE utf8_unicode_ci DEFAULT NULL,
  `cdna` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`acc_num`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hg19_coords`
--

DROP TABLE IF EXISTS `hg19_coords`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hg19_coords` (
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `chromosome` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `strand` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `coordSTART` int(11) DEFAULT NULL,
  `coordEND` int(11) DEFAULT NULL,
  `upstreamFLANK` varchar(30) COLLATE utf8_unicode_ci DEFAULT NULL,
  `downstreamFLANK` varchar(30) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hg38_coords`
--

DROP TABLE IF EXISTS `hg38_coords`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hg38_coords` (
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `chromosome` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `strand` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `coordSTART` int(11) DEFAULT NULL,
  `coordEND` int(11) DEFAULT NULL,
  `upstreamFLANK` varchar(30) COLLATE utf8_unicode_ci DEFAULT NULL,
  `downstreamFLANK` varchar(30) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hgmd_hg19_vcf`
--

DROP TABLE IF EXISTS `hgmd_hg19_vcf`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hgmd_hg19_vcf` (
  `chrom` varchar(45) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `pos` int(11) NOT NULL DEFAULT '0',
  `id` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `ref` varchar(255) COLLATE utf8_unicode_ci DEFAULT NULL,
  `alt` varchar(255) COLLATE utf8_unicode_ci DEFAULT NULL,
  `qual` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `filter` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `info` varchar(255) COLLATE utf8_unicode_ci DEFAULT NULL,
  `score` tinyint(3) unsigned DEFAULT NULL,
  PRIMARY KEY (`chrom`,`pos`,`id`),
  KEY `hgmd_score` (`score`),
  KEY `hgmd_id` (`id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hgmd_hg38_vcf`
--

DROP TABLE IF EXISTS `hgmd_hg38_vcf`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hgmd_hg38_vcf` (
  `chrom` varchar(2) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `pos` int(11) NOT NULL DEFAULT '0',
  `id` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `ref` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `alt` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `qual` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `filter` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `info` varchar(255) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`chrom`,`pos`,`id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `history`
--

DROP TABLE IF EXISTS `history`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `history` (
  `histID` int(11) NOT NULL DEFAULT '0',
  `acc_num` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `tname` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `colname` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `beforeUpd` text COLLATE utf8_unicode_ci,
  `afterUpd` text COLLATE utf8_unicode_ci,
  `updDate` date DEFAULT NULL,
  PRIMARY KEY (`histID`),
  KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `indel`
--

DROP TABLE IF EXISTS `indel`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `indel` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `wildtype` varchar(65) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `insertion` varchar(25) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `codon` int(11) NOT NULL DEFAULT '0',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`wildtype`,`insertion`,`codon`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `insertion`
--

DROP TABLE IF EXISTS `insertion`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `insertion` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `insertion` varchar(65) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `codon` int(11) NOT NULL DEFAULT '0',
  `nucleotide` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`insertion`,`codon`,`nucleotide`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `journal`
--

DROP TABLE IF EXISTS `journal`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `journal` (
  `journal` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `allname` varchar(200) COLLATE utf8_unicode_ci DEFAULT NULL,
  `url` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`journal`),
  KEY `fullname` (`fullname`),
  KEY `allname` (`allname`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `links`
--

DROP TABLE IF EXISTS `links`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `links` (
  `xorder` varchar(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `url` varchar(150) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`gene`,`url`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `markname`
--

DROP TABLE IF EXISTS `markname`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `markname` (
  `gene_id` int(11) NOT NULL DEFAULT '0',
  `genesym` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `chrom` varchar(25) COLLATE utf8_unicode_ci DEFAULT NULL,
  `genename` varchar(200) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gdbid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `omimid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `entrezID` int(11) DEFAULT NULL,
  `unigeneID` int(11) DEFAULT NULL,
  `hgncID` int(11) DEFAULT NULL,
  `hprdID` smallint(5) unsigned zerofill DEFAULT NULL,
  `swiss_acc` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene_date` date DEFAULT NULL,
  `svar` int(11) DEFAULT NULL,
  `gadb` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `findbase` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `cosmic` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `atlas` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `instruct` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mupit` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`gene_id`),
  UNIQUE KEY `entrezID` (`entrezID`),
  UNIQUE KEY `unigeneID` (`unigeneID`),
  UNIQUE KEY `hgncID` (`hgncID`),
  UNIQUE KEY `hprdID` (`hprdID`),
  UNIQUE KEY `genesym` (`genesym`),
  KEY `svar` (`svar`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutall`
--

DROP TABLE IF EXISTS `mutall`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutall` (
  `acc_num` varchar(10) DEFAULT NULL,
  `disease` varchar(65) DEFAULT NULL,
  `tag` enum('DP','FP','DFP','DM','FTV') DEFAULT NULL,
  KEY `acc_num_mutall` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutation`
--

DROP TABLE IF EXISTS `mutation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutation` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `base` varchar(9) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `amino` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `codon` int(11) NOT NULL DEFAULT '0',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`codon`,`base`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mutnomen`
--

DROP TABLE IF EXISTS `mutnomen`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mutnomen` (
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `refCORE` varchar(14) COLLATE utf8_unicode_ci DEFAULT NULL,
  `refVER` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `cSTART` varchar(11) COLLATE utf8_unicode_ci DEFAULT NULL,
  `ivsSTART` int(11) DEFAULT NULL,
  `cEND` varchar(11) COLLATE utf8_unicode_ci DEFAULT NULL,
  `ivsEND` int(11) DEFAULT NULL,
  `wildBASE` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutBASE` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `protCORE` varchar(14) COLLATE utf8_unicode_ci DEFAULT NULL,
  `protVER` varchar(2) COLLATE utf8_unicode_ci DEFAULT NULL,
  `wildAMINO` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `mutAMINO` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `codon` int(11) DEFAULT NULL,
  `hgvs` varchar(60) COLLATE utf8_unicode_ci DEFAULT NULL,
  `hgvsall` varchar(120) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`acc_num`),
  FULLTEXT KEY `HGVSALL` (`hgvsall`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `prom`
--

DROP TABLE IF EXISTS `prom`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `prom` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `base` varchar(65) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `location` int(11) NOT NULL DEFAULT '0',
  `locref` varchar(5) COLLATE utf8_unicode_ci DEFAULT NULL,
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`base`,`location`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `refseqgene`
--

DROP TABLE IF EXISTS `refseqgene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `refseqgene` (
  `entrez_id` int(11) NOT NULL DEFAULT '0',
  `ngseq` varchar(50) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`entrez_id`,`ngseq`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sec`
--

DROP TABLE IF EXISTS `sec`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sec` (
  `code3` char(4) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `helix` int(11) DEFAULT NULL,
  `sheet` int(11) DEFAULT NULL,
  `turn` int(11) DEFAULT NULL,
  `helixD` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `sheetD` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `p1` decimal(4,3) DEFAULT NULL,
  `p2` decimal(4,3) DEFAULT NULL,
  `p3` decimal(4,3) DEFAULT NULL,
  `p4` decimal(4,3) DEFAULT NULL,
  PRIMARY KEY (`code3`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `splice`
--

DROP TABLE IF EXISTS `splice`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `splice` (
  `disease` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  `gene` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `ivs` varchar(4) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `type` char(2) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `base` char(3) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `location` int(11) NOT NULL DEFAULT '0',
  `tag` enum('DP','FP','DFP','DM','DM?','FTV','R') COLLATE utf8_unicode_ci DEFAULT NULL,
  `author` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `journal` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `fullname` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `vol` varchar(6) COLLATE utf8_unicode_ci DEFAULT NULL,
  `page` varchar(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `year` year(4) DEFAULT NULL,
  `pmid` varchar(8) COLLATE utf8_unicode_ci DEFAULT NULL,
  `comments` text COLLATE utf8_unicode_ci,
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `new_date` date DEFAULT NULL,
  PRIMARY KEY (`gene`,`ivs`,`type`,`base`,`location`),
  UNIQUE KEY `acc_num` (`acc_num`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2019-11-25 13:57:14
