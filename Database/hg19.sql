-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: hg19
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
) ENGINE=MyISAM AUTO_INCREMENT=40956387 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
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
  KEY `cadd_tmp_key` (`chrom`,`start`,`ref`,`alt`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ccdsKgMap`
--

DROP TABLE IF EXISTS `ccdsKgMap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ccdsKgMap` (
  `ccdsId` varchar(32) NOT NULL,
  `geneId` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `cdsSimilarity` float NOT NULL,
  KEY `ccdsId` (`ccdsId`),
  KEY `geneId` (`geneId`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ccds_2007`
--

DROP TABLE IF EXISTS `ccds_2007`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ccds_2007` (
  `idccds` int(11) NOT NULL AUTO_INCREMENT,
  `chromosome` varchar(10) DEFAULT NULL,
  `g_accession` varchar(20) DEFAULT NULL,
  `gene` varchar(20) DEFAULT NULL,
  `gene_id` int(11) DEFAULT NULL,
  `ccds_id` varchar(30) DEFAULT NULL,
  `ccds_status` varchar(40) DEFAULT NULL,
  `cds_strand` char(1) DEFAULT NULL,
  `cds_from` int(11) DEFAULT NULL,
  `cds_to` int(11) DEFAULT NULL,
  `cds_locations` longblob,
  PRIMARY KEY (`idccds`)
) ENGINE=InnoDB AUTO_INCREMENT=36119 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ccds_current`
--

DROP TABLE IF EXISTS `ccds_current`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ccds_current` (
  `idccds` int(11) NOT NULL AUTO_INCREMENT,
  `chromosome` varchar(10) DEFAULT NULL,
  `g_accession` varchar(20) DEFAULT NULL,
  `gene` varchar(20) DEFAULT NULL,
  `gene_id` int(11) DEFAULT NULL,
  `ccds_id` varchar(30) DEFAULT NULL,
  `ccds_status` varchar(40) DEFAULT NULL,
  `cds_strand` char(1) DEFAULT NULL,
  `cds_from` int(11) DEFAULT NULL,
  `cds_to` int(11) DEFAULT NULL,
  `cds_locations` longblob,
  PRIMARY KEY (`idccds`)
) ENGINE=InnoDB AUTO_INCREMENT=27817 DEFAULT CHARSET=latin1;
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
  `path` varchar(255) NOT NULL,
  KEY `chromstart` (`chrom`,`start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dgv`
--

DROP TABLE IF EXISTS `dgv`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dgv` (
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`(16),`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
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
-- Table structure for table `ensGtp`
--

DROP TABLE IF EXISTS `ensGtp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ensGtp` (
  `gene` char(20) NOT NULL,
  `transcript` char(20) NOT NULL,
  `protein` char(24) NOT NULL,
  UNIQUE KEY `transcript` (`transcript`(19)),
  KEY `gene` (`gene`(19)),
  KEY `protein` (`protein`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ensemblToGeneName`
--

DROP TABLE IF EXISTS `ensemblToGeneName`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ensemblToGeneName` (
  `name` varchar(255) NOT NULL,
  `value` varchar(255) NOT NULL,
  PRIMARY KEY (`name`(15)),
  KEY `value` (`value`(22))
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
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
  `start` int(11) DEFAULT NULL,
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL,
  `filter` set('PASS','VQSR') DEFAULT NULL,
  `ea_homref` int(11) unsigned DEFAULT NULL,
  `ea_het` int(11) unsigned DEFAULT NULL,
  `ea_homalt` int(11) unsigned DEFAULT NULL,
  `aa_homref` int(11) unsigned DEFAULT NULL,
  `aa_het` int(11) unsigned DEFAULT NULL,
  `aa_homalt` int(11) unsigned DEFAULT NULL,
  `popmax_af` float DEFAULT NULL,
  PRIMARY KEY (`idexac`),
  UNIQUE KEY `uniqueexac` (`chrom`,`start`,`refallele`,`allele`)
) ENGINE=MyISAM AUTO_INCREMENT=275581146 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `evsOLD`
--

DROP TABLE IF EXISTS `evsOLD`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `evsOLD` (
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
) ENGINE=MyISAM AUTO_INCREMENT=10195873 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `evs_gnomad2`
--

DROP TABLE IF EXISTS `evs_gnomad2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `evs_gnomad2` (
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
) ENGINE=MyISAM AUTO_INCREMENT=260657407 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `evs_old_20180313`
--

DROP TABLE IF EXISTS `evs_old_20180313`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `evs_old_20180313` (
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
) ENGINE=MyISAM AUTO_INCREMENT=274056991 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
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
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`(16),`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gencodeV20chrM`
--

DROP TABLE IF EXISTS `gencodeV20chrM`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gencodeV20chrM` (
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL,
  `cdsEnd` int(10) unsigned NOT NULL,
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `score` int(11) DEFAULT NULL,
  `geneSymbol` varchar(255) NOT NULL,
  `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `exonFrames` longblob NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`),
  KEY `geneSymbol` (`geneSymbol`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gencodeV20chrM_cds`
--

DROP TABLE IF EXISTS `gencodeV20chrM_cds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gencodeV20chrM_cds` (
  `name` varchar(255) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genomeinabottleregions`
--

DROP TABLE IF EXISTS `genomeinabottleregions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomeinabottleregions` (
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`,`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genomeinabottleregions_geo`
--

DROP TABLE IF EXISTS `genomeinabottleregions_geo`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomeinabottleregions_geo` (
  `region` linestring NOT NULL,
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`,`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`),
  SPATIAL KEY `region` (`region`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genomicsuperdups`
--

DROP TABLE IF EXISTS `genomicsuperdups`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genomicsuperdups` (
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`(16),`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gnomadExomes`
--

DROP TABLE IF EXISTS `gnomadExomes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gnomadExomes` (
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
) ENGINE=MyISAM AUTO_INCREMENT=17002785 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `homopolymerregions`
--

DROP TABLE IF EXISTS `homopolymerregions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `homopolymerregions` (
  `chrom` varchar(45) NOT NULL,
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  `nucleotide` varchar(1) NOT NULL,
  KEY `chrom` (`chrom`,`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hpoGene`
--

DROP TABLE IF EXISTS `hpoGene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hpoGene` (
  `entrezgeneid` int(10) NOT NULL,
  `genesymbol` varchar(150) NOT NULL,
  `hpoTermName` varchar(150) NOT NULL,
  `hpoID` varchar(31) NOT NULL,
  KEY `hpogenesymbol` (`genesymbol`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `idmapping`
--

DROP TABLE IF EXISTS `idmapping`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `idmapping` (
  `acc` varchar(40) NOT NULL,
  `database` varchar(40) DEFAULT NULL,
  `externalid` varchar(40) DEFAULT NULL,
  KEY `idmapping_acc` (`acc`),
  KEY `idmapping_id` (`externalid`)
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
) ENGINE=MyISAM AUTO_INCREMENT=210644787 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `keggClass`
--

DROP TABLE IF EXISTS `keggClass`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `keggClass` (
  `classId` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `isA` int(11) DEFAULT NULL,
  PRIMARY KEY (`classId`),
  KEY `keggClassisA` (`isA`),
  CONSTRAINT `keggClassisA` FOREIGN KEY (`isA`) REFERENCES `keggClass` (`classId`)
) ENGINE=InnoDB AUTO_INCREMENT=276 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `keggOntology`
--

DROP TABLE IF EXISTS `keggOntology`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `keggOntology` (
  `koId` varchar(40) NOT NULL,
  `definition` varchar(255) DEFAULT NULL,
  `enzymeId` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`koId`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `keggOntology2keggClass`
--

DROP TABLE IF EXISTS `keggOntology2keggClass`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `keggOntology2keggClass` (
  `koId` varchar(40) NOT NULL DEFAULT '',
  `classId` int(11) NOT NULL DEFAULT '0',
  `ref` varchar(40) DEFAULT NULL COMMENT 'The actual reference of the class in the ontology. Useful when its of type PATH.',
  PRIMARY KEY (`koId`,`classId`),
  KEY `keggOntology2keggClass_keggClass` (`classId`),
  CONSTRAINT `keggOntology2keggClass_keggClass` FOREIGN KEY (`classId`) REFERENCES `keggClass` (`classId`),
  CONSTRAINT `keggOntology2keggClass_keggOntology` FOREIGN KEY (`koId`) REFERENCES `keggOntology` (`koId`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `keggPathways`
--

DROP TABLE IF EXISTS `keggPathways`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `keggPathways` (
  `keggId` varchar(20) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`keggId`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
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
-- Table structure for table `kgXref2keggOntology`
--

DROP TABLE IF EXISTS `kgXref2keggOntology`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `kgXref2keggOntology` (
  `geneId` varchar(40) NOT NULL DEFAULT '',
  `koId` varchar(40) NOT NULL DEFAULT '',
  PRIMARY KEY (`geneId`,`koId`),
  KEY `kgXref2keggOntology_keggOntology` (`koId`),
  CONSTRAINT `kgXref2keggOntology_keggOntology` FOREIGN KEY (`koId`) REFERENCES `keggOntology` (`koId`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `kgXref_old`
--

DROP TABLE IF EXISTS `kgXref_old`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `kgXref_old` (
  `kgID` varchar(40) NOT NULL,
  `mRNA` varchar(40) DEFAULT NULL,
  `spID` varchar(40) DEFAULT NULL,
  `spDisplayID` varchar(40) DEFAULT NULL,
  `geneSymbol` varchar(40) DEFAULT NULL,
  `refseq` varchar(40) DEFAULT NULL,
  `protAcc` varchar(40) DEFAULT NULL,
  `description` longblob NOT NULL,
  `keggID` varchar(40) DEFAULT NULL,
  KEY `kgID` (`kgID`),
  KEY `mRNA` (`mRNA`),
  KEY `spID` (`spID`),
  KEY `spDisplayID` (`spDisplayID`),
  KEY `geneSymbol` (`geneSymbol`),
  KEY `refseq` (`refseq`),
  KEY `protAcc` (`protAcc`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownCanonical`
--

DROP TABLE IF EXISTS `knownCanonical`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownCanonical` (
  `chrom` varchar(255) NOT NULL DEFAULT '',
  `chromStart` int(11) NOT NULL DEFAULT '0',
  `chromEnd` int(11) NOT NULL DEFAULT '0',
  `clusterId` int(11) NOT NULL DEFAULT '0',
  `transcript` varchar(255) NOT NULL DEFAULT '',
  `protein` varchar(255) NOT NULL DEFAULT '',
  UNIQUE KEY `clusterId` (`clusterId`),
  KEY `transcript` (`transcript`(12)),
  KEY `protein` (`protein`(12))
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
  KEY `align` (`alignID`),
  KEY `knowngene_chromend` (`chrom`,`txEnd`),
  KEY `knowngene_chrompos` (`chrom`,`txStart`,`txEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownGenePep`
--

DROP TABLE IF EXISTS `knownGenePep`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGenePep` (
  `name` varchar(255) NOT NULL DEFAULT '',
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
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownGene_cds_old`
--

DROP TABLE IF EXISTS `knownGene_cds_old`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGene_cds_old` (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `knownGene_cds_test`
--

DROP TABLE IF EXISTS `knownGene_cds_test`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `knownGene_cds_test` (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lincRNAsTranscripts`
--

DROP TABLE IF EXISTS `lincRNAsTranscripts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lincRNAsTranscripts` (
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL,
  `cdsEnd` int(10) unsigned NOT NULL,
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`),
  KEY `lincrna_chrompos` (`chrom`,`txStart`,`txEnd`),
  KEY `lincrna_chromend` (`chrom`,`txEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary table structure for view `lincRNAsTranscriptsSymbol`
--

DROP TABLE IF EXISTS `lincRNAsTranscriptsSymbol`;
/*!50001 DROP VIEW IF EXISTS `lincRNAsTranscriptsSymbol`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `lincRNAsTranscriptsSymbol` (
  `bin` tinyint NOT NULL,
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
  `geneSymbol` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `lincRNAsTranscripts_cds`
--

DROP TABLE IF EXISTS `lincRNAsTranscripts_cds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lincRNAsTranscripts_cds` (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Temporary table structure for view `longest_gene`
--

DROP TABLE IF EXISTS `longest_gene`;
/*!50001 DROP VIEW IF EXISTS `longest_gene`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `longest_gene` (
  `geneSymbol` tinyint NOT NULL,
  `longest_transcript` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `longest_peptide`
--

DROP TABLE IF EXISTS `longest_peptide`;
/*!50001 DROP VIEW IF EXISTS `longest_peptide`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `longest_peptide` (
  `geneSymbol` tinyint NOT NULL,
  `longest_transcript` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `longest_transcript`
--

DROP TABLE IF EXISTS `longest_transcript`;
/*!50001 DROP VIEW IF EXISTS `longest_transcript`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `longest_transcript` (
  `geneSymbol` tinyint NOT NULL,
  `longest_transcript` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Table structure for table `lowcomplexityregions`
--

DROP TABLE IF EXISTS `lowcomplexityregions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `lowcomplexityregions` (
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`,`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `miRNA`
--

DROP TABLE IF EXISTS `miRNA`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `miRNA` (
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL,
  `cdsEnd` int(10) unsigned NOT NULL,
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `geneSymbol` varchar(255) DEFAULT NULL,
  `type` varchar(255) NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`),
  KEY `mirna_chrompos` (`chrom`,`txStart`,`txEnd`),
  KEY `mirna_chromend` (`chrom`,`txEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `miRNA_cds`
--

DROP TABLE IF EXISTS `miRNA_cds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `miRNA_cds` (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mitodomains`
--

DROP TABLE IF EXISTS `mitodomains`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mitodomains` (
  `idmitodomains` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL DEFAULT '',
  `chrom` varchar(45) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  PRIMARY KEY (`idmitodomains`),
  KEY `name` (`name`),
  KEY `chrom` (`chrom`(16),`start`)
) ENGINE=InnoDB AUTO_INCREMENT=18771 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mitogb`
--

DROP TABLE IF EXISTS `mitogb`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mitogb` (
  `idgenbank` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL,
  `haplogroup` varchar(255) NOT NULL DEFAULT '',
  `haplogroup_split` varchar(255) NOT NULL DEFAULT '',
  `af` float NOT NULL DEFAULT '0',
  PRIMARY KEY (`idgenbank`),
  KEY `chromstart` (`chrom`,`start`),
  KEY `chromend` (`chrom`,`end`),
  KEY `af` (`af`)
) ENGINE=InnoDB AUTO_INCREMENT=77252 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mitomap`
--

DROP TABLE IF EXISTS `mitomap`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mitomap` (
  `chrom` varchar(45) NOT NULL DEFAULT '',
  `start` int(11) NOT NULL DEFAULT '0',
  `end` int(11) unsigned NOT NULL DEFAULT '0',
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL DEFAULT '',
  `af` float DEFAULT '0',
  `reference_links` varchar(255) DEFAULT NULL,
  `disease` varchar(255) DEFAULT '',
  `disease_homoplasmy` varchar(1) DEFAULT NULL,
  `disease_heteroplasmy` varchar(1) DEFAULT NULL,
  `disease_status` varchar(45) DEFAULT NULL,
  `disease_reference_links` varchar(255) DEFAULT NULL,
  `tissue` varchar(255) DEFAULT '',
  `tissue_homoplasmy` varchar(1) DEFAULT NULL,
  `tissue_heteroplasmy` varchar(1) DEFAULT NULL,
  `tissue_reference_links` varchar(255) DEFAULT NULL,
  `idmitomap` int(11) unsigned NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`idmitomap`),
  UNIQUE KEY `uniquemitomap` (`chrom`,`start`,`end`,`refallele`,`allele`),
  KEY `chromstart` (`chrom`,`start`),
  KEY `chromend` (`chrom`,`end`)
) ENGINE=InnoDB AUTO_INCREMENT=11499 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `phylopgenes`
--

DROP TABLE IF EXISTS `phylopgenes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phylopgenes` (
  `idphylopgenes` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `genesymbol` varchar(150) NOT NULL,
  `phylopscore` float NOT NULL,
  PRIMARY KEY (`idphylopgenes`),
  KEY `higenesymbol` (`genesymbol`)
) ENGINE=InnoDB AUTO_INCREMENT=19593 DEFAULT CHARSET=latin1;
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
-- Table structure for table `proteinaliases`
--

DROP TABLE IF EXISTS `proteinaliases`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `proteinaliases` (
  `string_protein_id` char(20) NOT NULL DEFAULT '',
  `alias` varchar(255) NOT NULL DEFAULT '',
  `source` varchar(255) NOT NULL DEFAULT '',
  KEY `string_protein_id` (`string_protein_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `proteinlinks`
--

DROP TABLE IF EXISTS `proteinlinks`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `proteinlinks` (
  `protein1` char(20) NOT NULL DEFAULT '',
  `protein2` char(20) NOT NULL DEFAULT '',
  `score` int(11) NOT NULL DEFAULT '0',
  KEY `protein1` (`protein1`),
  KEY `protein2` (`protein2`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `refGene`
--

DROP TABLE IF EXISTS `refGene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `refGene` (
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL,
  `cdsEnd` int(10) unsigned NOT NULL,
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `score` int(11) DEFAULT NULL,
  `geneSymbol` varchar(255) NOT NULL,
  `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `exonFrames` longblob NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`),
  KEY `geneSymbol` (`geneSymbol`),
  KEY `chrompos` (`chrom`,`txStart`,`txEnd`),
  KEY `chromend` (`chrom`,`txEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `refGene_cds`
--

DROP TABLE IF EXISTS `refGene_cds`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `refGene_cds` (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `refGene_cds_old`
--

DROP TABLE IF EXISTS `refGene_cds_old`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `refGene_cds_old` (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `refGene_old`
--

DROP TABLE IF EXISTS `refGene_old`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `refGene_old` (
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
  `geneSymbol` varchar(40) DEFAULT NULL,
  KEY `name` (`name`),
  KEY `chrom` (`chrom`(16),`txStart`),
  KEY `chrom_2` (`chrom`(16),`txEnd`),
  KEY `protein` (`proteinID`(16)),
  KEY `align` (`alignID`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `rmsk`
--

DROP TABLE IF EXISTS `rmsk`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `rmsk` (
  `bin` smallint(5) unsigned NOT NULL DEFAULT '0',
  `swScore` int(10) unsigned NOT NULL DEFAULT '0',
  `milliDiv` int(10) unsigned NOT NULL DEFAULT '0',
  `milliDel` int(10) unsigned NOT NULL DEFAULT '0',
  `milliIns` int(10) unsigned NOT NULL DEFAULT '0',
  `genoName` varchar(255) NOT NULL DEFAULT '',
  `genoStart` int(10) unsigned NOT NULL DEFAULT '0',
  `genoEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `genoLeft` int(11) NOT NULL DEFAULT '0',
  `strand` char(1) NOT NULL DEFAULT '',
  `repName` varchar(255) NOT NULL DEFAULT '',
  `repClass` varchar(255) NOT NULL DEFAULT '',
  `repFamily` varchar(255) NOT NULL DEFAULT '',
  `repStart` int(11) NOT NULL DEFAULT '0',
  `repEnd` int(11) NOT NULL DEFAULT '0',
  `repLeft` int(11) NOT NULL DEFAULT '0',
  `id` char(1) NOT NULL DEFAULT '',
  KEY `genoName` (`genoName`(14),`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sec_ac`
--

DROP TABLE IF EXISTS `sec_ac`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sec_ac` (
  `secondary` varchar(6) NOT NULL DEFAULT '',
  `primary` varchar(6) NOT NULL DEFAULT '',
  KEY `sec` (`secondary`),
  KEY `pri` (`primary`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
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
-- Table structure for table `snp130`
--

DROP TABLE IF EXISTS `snp130`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp130` (
  `bin` smallint(5) unsigned NOT NULL DEFAULT '0',
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `chromStart` int(10) unsigned NOT NULL DEFAULT '0',
  `chromEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `name` varchar(15) NOT NULL DEFAULT '',
  `score` smallint(5) unsigned NOT NULL DEFAULT '0',
  `strand` enum('+','-') DEFAULT NULL,
  `refNCBI` blob NOT NULL,
  `refUCSC` blob NOT NULL,
  `observed` varchar(255) NOT NULL DEFAULT '',
  `molType` enum('unknown','genomic','cDNA') DEFAULT NULL,
  `class` enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion') NOT NULL DEFAULT 'unknown',
  `valid` set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes') NOT NULL DEFAULT 'unknown',
  `avHet` float NOT NULL DEFAULT '0',
  `avHetSE` float NOT NULL DEFAULT '0',
  `func` set('unknown','coding-synon','intron','near-gene-3','near-gene-5','nonsense','missense','frameshift','untranslated-3','untranslated-5') NOT NULL DEFAULT 'unknown',
  `locType` enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion') DEFAULT NULL,
  `weight` int(10) unsigned NOT NULL DEFAULT '0',
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snp131`
--

DROP TABLE IF EXISTS `snp131`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp131` (
  `bin` smallint(5) unsigned NOT NULL DEFAULT '0',
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `chromStart` int(10) unsigned NOT NULL DEFAULT '0',
  `chromEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `name` varchar(15) NOT NULL DEFAULT '',
  `score` smallint(5) unsigned NOT NULL DEFAULT '0',
  `strand` enum('+','-') DEFAULT NULL,
  `refNCBI` blob NOT NULL,
  `refUCSC` blob NOT NULL,
  `observed` varchar(255) NOT NULL DEFAULT '',
  `molType` enum('unknown','genomic','cDNA') DEFAULT NULL,
  `class` enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion') NOT NULL DEFAULT 'unknown',
  `valid` set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes') NOT NULL DEFAULT 'unknown',
  `avHet` float NOT NULL DEFAULT '0',
  `avHetSE` float NOT NULL DEFAULT '0',
  `func` set('unknown','coding-synon','intron','coding-synonymy-unknown','near-gene-3','near-gene-5','nonsense','missense','frameshift','cds-indel','untranslated-3','untranslated-5','splice-3','splice-5') NOT NULL DEFAULT 'unknown',
  `locType` enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion') DEFAULT NULL,
  `weight` int(10) unsigned NOT NULL DEFAULT '0',
  `clinical` varchar(25) DEFAULT '',
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`),
  KEY `chrom2` (`chrom`,`chromStart`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snp131OrthoPt2Pa2Rm2`
--

DROP TABLE IF EXISTS `snp131OrthoPt2Pa2Rm2`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp131OrthoPt2Pa2Rm2` (
  `bin` smallint(6) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `humanObserved` varchar(255) NOT NULL,
  `humanAllele` char(1) NOT NULL,
  `humanStrand` char(1) NOT NULL,
  `chimpChrom` varchar(255) NOT NULL,
  `chimpStart` int(10) unsigned NOT NULL,
  `chimpEnd` int(10) unsigned NOT NULL,
  `chimpAllele` char(1) NOT NULL,
  `chimpStrand` char(1) NOT NULL,
  `orangChrom` varchar(255) NOT NULL,
  `orangStart` int(10) unsigned NOT NULL,
  `orangEnd` int(10) unsigned NOT NULL,
  `orangAllele` char(1) NOT NULL,
  `orangStrand` char(1) NOT NULL,
  `macaqueChrom` varchar(255) NOT NULL,
  `macaqueStart` int(10) unsigned NOT NULL,
  `macaqueEnd` int(10) unsigned NOT NULL,
  `macaqueAllele` char(1) NOT NULL,
  `macaqueStrand` char(1) NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snp132`
--

DROP TABLE IF EXISTS `snp132`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp132` (
  `bin` smallint(5) unsigned NOT NULL DEFAULT '0',
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `chromStart` int(10) unsigned NOT NULL DEFAULT '0',
  `chromEnd` int(10) unsigned NOT NULL DEFAULT '0',
  `name` varchar(15) NOT NULL DEFAULT '',
  `score` smallint(5) unsigned NOT NULL DEFAULT '0',
  `strand` enum('+','-') DEFAULT NULL,
  `refNCBI` blob NOT NULL,
  `refUCSC` blob NOT NULL,
  `observed` varchar(255) NOT NULL DEFAULT '',
  `molType` enum('unknown','genomic','cDNA') DEFAULT NULL,
  `class` enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion') NOT NULL DEFAULT 'unknown',
  `valid` set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes') NOT NULL DEFAULT 'unknown',
  `avHet` float NOT NULL DEFAULT '0',
  `avHetSE` float NOT NULL DEFAULT '0',
  `func` set('unknown','coding-synon','intron','coding-synonymy-unknown','near-gene-3','near-gene-5','nonsense','missense','frameshift','cds-indel','untranslated-3','untranslated-5','splice-3','splice-5') NOT NULL DEFAULT 'unknown',
  `locType` enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion') DEFAULT NULL,
  `weight` int(10) unsigned NOT NULL DEFAULT '0',
  `clinical` varchar(25) DEFAULT '',
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`),
  KEY `chrom2` (`chrom`,`chromStart`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snp135`
--

DROP TABLE IF EXISTS `snp135`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snp135` (
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
  `molType` enum('unknown','genomic','cDNA') NOT NULL,
  `class` enum('unknown','single','in-del','het','microsatellite','named','mixed','mnp','insertion','deletion') NOT NULL,
  `valid` set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes') NOT NULL,
  `avHet` float NOT NULL,
  `avHetSE` float NOT NULL,
  `func` set('unknown','coding-synon','intron','near-gene-3','near-gene-5','nonsense','missense','stop-loss','frameshift','cds-indel','untranslated-3','untranslated-5','splice-5') NOT NULL,
  `locType` enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion') NOT NULL,
  `weight` int(10) unsigned NOT NULL,
  `exceptions` set('RefAlleleMismatch','RefAlleleRevComp','DuplicateObserved','MixedObserved','FlankMismatchGenomeLonger','FlankMismatchGenomeEqual','FlankMismatchGenomeShorter','NamedDeletionZeroSpan','NamedInsertionNonzeroSpan','SingleClassLongerSpan','SingleClassZeroSpan','SingleClassTriAllelic','SingleClassQuadAllelic','ObservedWrongFormat','ObservedTooLong','ObservedContainsIupac','ObservedMismatch','MultipleAlignments','NonIntegerChromCount','AlleleFreqSumNot1','SingleAlleleFreq','InconsistentAlleles') NOT NULL,
  `submitterCount` smallint(5) unsigned NOT NULL,
  `submitters` longblob NOT NULL,
  `alleleFreqCount` smallint(5) unsigned NOT NULL,
  `alleles` longblob NOT NULL,
  `alleleNs` longblob NOT NULL,
  `alleleFreqs` longblob NOT NULL,
  `bitfields` set('clinically-assoc','maf-5-some-pop','maf-5-all-pops','has-omim-omia','microattr-tpa','submitted-by-lsdb','genotype-conflict','rs-cluster-nonoverlapping-alleles','observed-mismatch') NOT NULL,
  `clinical` varchar(25) DEFAULT '',
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`),
  KEY `chrom2` (`chrom`,`chromStart`)
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
  `molType` enum('unknown','genomic','cDNA') NOT NULL,
  `class` enum('unknown','single','in-del','microsatellite','named','mnp','insertion','deletion') NOT NULL,
  `valid` set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes') NOT NULL,
  `avHet` float NOT NULL,
  `avHetSE` float NOT NULL,
  `func` set('unknown','coding-synon','intron','near-gene-3','near-gene-5','ncRNA','nonsense','missense','stop-loss','frameshift','cds-indel','untranslated-3','untranslated-5','splice-3','splice-5') NOT NULL,
  `locType` enum('range','exact','between','rangeInsertion','rangeSubstitution','rangeDeletion','fuzzy') NOT NULL,
  `weight` int(10) unsigned NOT NULL,
  `exceptions` set('RefAlleleMismatch','RefAlleleRevComp','DuplicateObserved','MixedObserved','FlankMismatchGenomeLonger','FlankMismatchGenomeEqual','FlankMismatchGenomeShorter','SingleClassLongerSpan','SingleClassZeroSpan','SingleClassTriAllelic','SingleClassQuadAllelic','ObservedWrongFormat','ObservedTooLong','ObservedContainsIupac','ObservedMismatch','MultipleAlignments','NonIntegerChromCount','AlleleFreqSumNot1','SingleAlleleFreq','InconsistentAlleles') NOT NULL,
  `submitterCount` smallint(5) unsigned NOT NULL,
  `submitters` longblob NOT NULL,
  `alleleFreqCount` smallint(5) unsigned NOT NULL,
  `alleles` longblob NOT NULL,
  `alleleNs` longblob NOT NULL,
  `alleleFreqs` longblob NOT NULL,
  `bitfields` set('clinically-assoc','maf-5-some-pop','maf-5-all-pops','has-omim-omia','microattr-tpa','submitted-by-lsdb','genotype-conflict','rs-cluster-nonoverlapping-alleles','observed-mismatch') NOT NULL,
  KEY `name` (`name`),
  KEY `chrom` (`chrom`,`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tgenomes`
--

DROP TABLE IF EXISTS `tgenomes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tgenomes` (
  `chrom` varchar(45) NOT NULL DEFAULT '',
  `pos` int(11) NOT NULL DEFAULT '0',
  `rssnp` varchar(15) DEFAULT '',
  `ref` varchar(255) DEFAULT '',
  `alt` varchar(255) DEFAULT '',
  `quality` varchar(15) DEFAULT '',
  `pass` varchar(15) DEFAULT '',
  `dp` int(10) unsigned NOT NULL DEFAULT '0',
  `af` float NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`,`pos`),
  KEY `chromref` (`chrom`,`pos`,`ref`,`alt`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tgenomes_old`
--

DROP TABLE IF EXISTS `tgenomes_old`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tgenomes_old` (
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `rssnp` varchar(15) DEFAULT '',
  `ref` varchar(15) DEFAULT '',
  `alt` varchar(15) DEFAULT '',
  `dummy` varchar(15) DEFAULT '',
  `pass` varchar(15) DEFAULT '',
  `dp` int(10) unsigned NOT NULL DEFAULT '0',
  `af` float NOT NULL DEFAULT '0',
  `cb` varchar(15) DEFAULT '',
  `eur` float DEFAULT NULL,
  KEY `chrom` (`chrom`,`pos`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `tgenomes_sv`
--

DROP TABLE IF EXISTS `tgenomes_sv`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `tgenomes_sv` (
  `chrom` varchar(31) NOT NULL DEFAULT '',
  `pos` int(10) unsigned NOT NULL DEFAULT '0',
  `end` int(10) unsigned NOT NULL DEFAULT '0',
  `svtype` enum('CNV','DEL','DEL_ALU','DEL_HERV','DEL_LINE1','DEL_SVA','DUP','INS','INV') NOT NULL,
  `af` float NOT NULL DEFAULT '0',
  KEY `chrom` (`chrom`,`pos`,`end`),
  KEY `chrom_1` (`chrom`(16),`pos`),
  KEY `chrom_2` (`chrom`(16),`end`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `wgEncodeAwgSegmentationCombinedH1hesc`
--

DROP TABLE IF EXISTS `wgEncodeAwgSegmentationCombinedH1hesc`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `wgEncodeAwgSegmentationCombinedH1hesc` (
  `bin` smallint(5) unsigned NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `chromStart` int(10) unsigned NOT NULL,
  `chromEnd` int(10) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `score` int(10) unsigned NOT NULL,
  `strand` char(1) NOT NULL,
  `thickStart` int(10) unsigned NOT NULL,
  `thickEnd` int(10) unsigned NOT NULL,
  `reserved` int(10) unsigned NOT NULL,
  KEY `name` (`name`(16)),
  KEY `chrom` (`chrom`(14),`bin`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `wgEncodeGencodeAttrsV12`
--

DROP TABLE IF EXISTS `wgEncodeGencodeAttrsV12`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `wgEncodeGencodeAttrsV12` (
  `geneId` varchar(255) NOT NULL,
  `geneName` varchar(255) NOT NULL,
  `geneType` varchar(255) NOT NULL,
  `geneStatus` varchar(255) NOT NULL,
  `transcriptId` varchar(255) NOT NULL,
  `transcriptName` varchar(255) NOT NULL,
  `transcriptType` varchar(255) NOT NULL,
  `transcriptStatus` varchar(255) NOT NULL,
  `havanaGeneId` varchar(255) NOT NULL,
  `havanaTranscriptId` varchar(255) NOT NULL,
  `ccdsId` varchar(255) NOT NULL,
  `level` int(11) NOT NULL,
  `transcriptClass` varchar(255) NOT NULL,
  PRIMARY KEY (`transcriptId`),
  KEY `geneId` (`geneId`),
  KEY `geneName` (`geneName`),
  KEY `havanaGeneId` (`havanaGeneId`),
  KEY `havanaTranscriptId` (`havanaTranscriptId`),
  KEY `ccdsId` (`ccdsId`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `wgEncodeGencodeBasicV12`
--

DROP TABLE IF EXISTS `wgEncodeGencodeBasicV12`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `wgEncodeGencodeBasicV12` (
  `bin` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `chrom` varchar(255) NOT NULL,
  `strand` char(1) NOT NULL,
  `txStart` int(10) unsigned NOT NULL,
  `txEnd` int(10) unsigned NOT NULL,
  `cdsStart` int(10) unsigned NOT NULL,
  `cdsEnd` int(10) unsigned NOT NULL,
  `exonCount` int(10) unsigned NOT NULL,
  `exonStarts` longblob NOT NULL,
  `exonEnds` longblob NOT NULL,
  `score` int(11) DEFAULT NULL,
  `name2` varchar(255) NOT NULL,
  `cdsStartStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `cdsEndStat` enum('none','unk','incmpl','cmpl') NOT NULL,
  `exonFrames` longblob NOT NULL,
  KEY `chrom` (`chrom`,`bin`),
  KEY `name` (`name`),
  KEY `name2` (`name2`)
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

--
-- Final view structure for view `lincRNAsTranscriptsSymbol`
--

/*!50001 DROP TABLE IF EXISTS `lincRNAsTranscriptsSymbol`*/;
/*!50001 DROP VIEW IF EXISTS `lincRNAsTranscriptsSymbol`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`wieland`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `lincRNAsTranscriptsSymbol` AS select `lincRNAsTranscripts`.`bin` AS `bin`,`lincRNAsTranscripts`.`name` AS `name`,`lincRNAsTranscripts`.`chrom` AS `chrom`,`lincRNAsTranscripts`.`strand` AS `strand`,`lincRNAsTranscripts`.`txStart` AS `txStart`,`lincRNAsTranscripts`.`txEnd` AS `txEnd`,`lincRNAsTranscripts`.`cdsStart` AS `cdsStart`,`lincRNAsTranscripts`.`cdsEnd` AS `cdsEnd`,`lincRNAsTranscripts`.`exonCount` AS `exonCount`,`lincRNAsTranscripts`.`exonStarts` AS `exonStarts`,`lincRNAsTranscripts`.`exonEnds` AS `exonEnds`,'' AS `geneSymbol` from `lincRNAsTranscripts` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `longest_gene`
--

/*!50001 DROP TABLE IF EXISTS `longest_gene`*/;
/*!50001 DROP VIEW IF EXISTS `longest_gene`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`wieland`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `longest_gene` AS select `longest_peptide`.`geneSymbol` AS `geneSymbol`,`longest_peptide`.`longest_transcript` AS `longest_transcript` from `longest_peptide` union select `longest_transcript`.`geneSymbol` AS `geneSymbol`,`longest_transcript`.`longest_transcript` AS `longest_transcript` from `longest_transcript` where (not(`longest_transcript`.`geneSymbol` in (select `longest_peptide`.`geneSymbol` from `longest_peptide`))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `longest_peptide`
--

/*!50001 DROP TABLE IF EXISTS `longest_peptide`*/;
/*!50001 DROP VIEW IF EXISTS `longest_peptide`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`wieland`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `longest_peptide` AS select distinct `g`.`geneSymbol` AS `geneSymbol`,(select `x`.`kgID` from (`kgXref` `x` join `knownGenePep` `p` on((`x`.`kgID` = `p`.`name`))) where (`g`.`geneSymbol` = `x`.`geneSymbol`) group by `g`.`geneSymbol` having max(length(`p`.`seq`))) AS `longest_transcript` from `kgXref` `g` group by `g`.`geneSymbol` having (`longest_transcript` is not null) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `longest_transcript`
--

/*!50001 DROP TABLE IF EXISTS `longest_transcript`*/;
/*!50001 DROP VIEW IF EXISTS `longest_transcript`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`wieland`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `longest_transcript` AS select distinct `g`.`geneSymbol` AS `geneSymbol`,(select `x`.`kgID` from (`kgXref` `x` join `knownGene` `p` on((`x`.`kgID` = `p`.`name`))) where (`g`.`geneSymbol` = `x`.`geneSymbol`) group by `g`.`geneSymbol` having max((`p`.`txEnd` - `p`.`txStart`))) AS `longest_transcript` from `kgXref` `g` group by `g`.`geneSymbol` having (`longest_transcript` is not null) */;
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
