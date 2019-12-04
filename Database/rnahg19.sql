-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: rnahg19
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
-- Table structure for table `additionalannotation`
--

DROP TABLE IF EXISTS `additionalannotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `additionalannotation` (
  `idadditionalannotation` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsnv` int(11) unsigned NOT NULL,
  `annotationname` enum('regulation') DEFAULT NULL,
  `idannotation` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idadditionalannotation`),
  KEY `additionalannotationidsnv` (`idsnv`),
  CONSTRAINT `additionalannotationidsnv` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=509587 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `decondition`
--

DROP TABLE IF EXISTS `decondition`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `decondition` (
  `idcondition` int(12) unsigned NOT NULL AUTO_INCREMENT,
  `clabel` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`idcondition`)
) ENGINE=MyISAM AUTO_INCREMENT=3 DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `deexperiment`
--

DROP TABLE IF EXISTS `deexperiment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `deexperiment` (
  `iddeexperiment` int(12) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(45) DEFAULT NULL,
  `tool` enum('deseq','deseq2','edger') NOT NULL,
  PRIMARY KEY (`iddeexperiment`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `deresult`
--

DROP TABLE IF EXISTS `deresult`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `deresult` (
  `idderesult` int(12) unsigned NOT NULL AUTO_INCREMENT,
  `idgene` int(11) unsigned NOT NULL,
  `iddeexperiment` int(12) unsigned NOT NULL,
  `log2FC` float DEFAULT NULL,
  `pvalue` float DEFAULT NULL,
  `padjusted` float DEFAULT NULL,
  PRIMARY KEY (`idderesult`),
  KEY `deresiddeexperiment` (`iddeexperiment`),
  KEY `deresidresult_index` (`idderesult`) USING BTREE,
  KEY `deresidgene_index` (`idgene`) USING BTREE
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `disease2gene`
--

DROP TABLE IF EXISTS `disease2gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `disease2gene` (
  `iddisease2gene` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `iddisease` int(11) unsigned NOT NULL,
  `idgene` int(11) unsigned NOT NULL,
  `class` int(11) DEFAULT NULL,
  PRIMARY KEY (`iddisease2gene`),
  UNIQUE KEY `iddisease2gene_UNIQUE` (`iddisease2gene`),
  UNIQUE KEY `diseasegene` (`iddisease`,`idgene`),
  KEY `d2giddisease` (`iddisease`),
  KEY `d2gidgene` (`idgene`),
  CONSTRAINT `d2giddisease` FOREIGN KEY (`iddisease`) REFERENCES `exomehg19`.`disease` (`iddisease`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `d2gidgene` FOREIGN KEY (`idgene`) REFERENCES `gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `exomestat`
--

DROP TABLE IF EXISTS `exomestat`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exomestat` (
  `idexomestat` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsample` int(11) unsigned NOT NULL,
  `duplicates` decimal(5,2) unsigned DEFAULT NULL,
  `reads` int(11) unsigned DEFAULT NULL,
  `mapped` int(11) unsigned DEFAULT NULL,
  `percentm` decimal(5,2) unsigned DEFAULT NULL,
  `seq` decimal(18,3) unsigned DEFAULT NULL,
  `onbait` decimal(5,2) unsigned DEFAULT NULL,
  `avgcov` decimal(7,2) unsigned DEFAULT NULL,
  `uncovered` decimal(5,2) unsigned DEFAULT NULL,
  `cov1x` decimal(5,2) unsigned DEFAULT NULL,
  `cov4x` decimal(5,2) unsigned DEFAULT NULL,
  `cov8x` decimal(5,2) unsigned DEFAULT NULL,
  `cov20x` decimal(5,2) unsigned DEFAULT NULL,
  `tstv` decimal(5,2) DEFAULT NULL,
  `idlibtype` int(11) unsigned NOT NULL,
  `sry` int(11) unsigned DEFAULT NULL,
  `idlibpair` int(11) unsigned NOT NULL,
  `mix` float DEFAULT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `avgcovstd` decimal(5,2) unsigned DEFAULT NULL,
  `mediancov` decimal(5,2) unsigned DEFAULT NULL,
  `mediancovstd` decimal(5,2) unsigned DEFAULT NULL,
  `exomedepthrsd` float DEFAULT NULL,
  `idassay` int(11) unsigned DEFAULT NULL,
  `mismatchrate` float DEFAULT NULL,
  `avgqual` float DEFAULT NULL,
  `avgquallast5` float DEFAULT NULL,
  `libcomplexity` float DEFAULT NULL,
  `q30fraction` float DEFAULT NULL,
  `opticalduplicates` decimal(5,2) unsigned DEFAULT NULL,
  `avgdiffdepth` decimal(5,2) DEFAULT NULL,
  `properlyp` decimal(5,2) unsigned DEFAULT NULL,
  PRIMARY KEY (`idexomestat`),
  UNIQUE KEY `exomeidsample` (`idlibtype`,`idlibpair`,`idsample`),
  KEY `exomestatsample` (`idsample`),
  KEY `exomestatlibpair` (`idlibpair`)
) ENGINE=InnoDB AUTO_INCREMENT=38882 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `funcrank`
--

DROP TABLE IF EXISTS `funcrank`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `funcrank` (
  `idfuncrank` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `func` varchar(50) DEFAULT NULL,
  `rank` smallint(5) unsigned DEFAULT NULL,
  PRIMARY KEY (`idfuncrank`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene`
--

DROP TABLE IF EXISTS `gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene` (
  `idgene` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `genesymbol` varchar(150) NOT NULL DEFAULT '',
  `omim` int(11) unsigned DEFAULT NULL,
  `blacklist` tinyint(3) unsigned DEFAULT '0',
  `nonsynpergene` int(10) unsigned NOT NULL,
  `maxpeplength` int(10) unsigned DEFAULT NULL,
  `meanreaddepth` float DEFAULT NULL,
  `readdepth20` float DEFAULT NULL,
  `delpergene` int(10) unsigned NOT NULL,
  `cgd` set('','AD','AR','XL','BG','Digenic','Methylation','Maternal','Paternal','Multigenic','PAR') DEFAULT NULL,
  PRIMARY KEY (`idgene`),
  UNIQUE KEY `genesymbol` (`genesymbol`),
  KEY `omim` (`omim`)
) ENGINE=InnoDB AUTO_INCREMENT=3563767 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `genebased`
--

DROP TABLE IF EXISTS `genebased`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `genebased` (
  `idgenebased` int(12) unsigned NOT NULL AUTO_INCREMENT,
  `idsample` int(11) unsigned NOT NULL,
  `idgene` int(3) unsigned NOT NULL,
  `readcount` int(10) unsigned DEFAULT NULL,
  `fpkm` double unsigned DEFAULT NULL,
  `date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`idgenebased`),
  UNIQUE KEY `idgeneidsample` (`idgene`,`idsample`),
  UNIQUE KEY `idsampleidgene` (`idsample`,`idgene`),
  KEY `idsample_index` (`idsample`) USING BTREE,
  KEY `idgene_index` (`idgene`) USING BTREE,
  CONSTRAINT `idgene` FOREIGN KEY (`idgene`) REFERENCES `exomevcf`.`gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `idsample` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=148102024 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `homozygosity`
--

DROP TABLE IF EXISTS `homozygosity`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `homozygosity` (
  `idhomozygosity` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsample` int(11) unsigned NOT NULL,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `count` int(11) unsigned NOT NULL,
  `score` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idhomozygosity`),
  UNIQUE KEY `homoidsample` (`idsample`,`chrom`,`start`),
  KEY `idsamplehomo` (`idsample`),
  KEY `chromstarthomo` (`chrom`,`start`),
  CONSTRAINT `idsamplehomo` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sample2deexperiment`
--

DROP TABLE IF EXISTS `sample2deexperiment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sample2deexperiment` (
  `idsample2deexperiment` int(12) unsigned NOT NULL AUTO_INCREMENT,
  `idsample` int(11) unsigned NOT NULL,
  `iddecondition` int(12) unsigned NOT NULL,
  `iddeexperiment` int(12) unsigned NOT NULL,
  PRIMARY KEY (`idsample2deexperiment`),
  KEY `s2deidsample` (`idsample`),
  KEY `s2deiddecondition` (`iddecondition`),
  KEY `s2deiddeexperiment` (`iddeexperiment`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 PAGE_CHECKSUM=1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snpeff`
--

DROP TABLE IF EXISTS `snpeff`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snpeff` (
  `idsnpeff` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsnv` int(11) unsigned NOT NULL,
  `idgene` int(11) unsigned DEFAULT NULL,
  `idsnpeffeffect` int(11) unsigned NOT NULL,
  `transcript` varchar(20) DEFAULT NULL,
  `effect_impact` enum('HIGH','MODERATE','LOW','MODIFIER') DEFAULT NULL,
  `functional_class` enum('NONE','SILENT','MISSENSE','NONSENSE') DEFAULT NULL,
  `codon_change` varchar(10) DEFAULT NULL,
  `amino_acid_change` varchar(100) DEFAULT NULL,
  `amino_acid_length` int(11) unsigned DEFAULT NULL,
  `exon_rank` int(11) unsigned DEFAULT NULL,
  PRIMARY KEY (`idsnpeff`),
  UNIQUE KEY `uniquesnpeff` (`idsnv`,`transcript`),
  KEY `snpeff_gene` (`idgene`),
  KEY `snpeff_snpeffeffect` (`idsnpeffeffect`),
  CONSTRAINT `snpeff_gene` FOREIGN KEY (`idgene`) REFERENCES `gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `snpeff_snpeffeffect` FOREIGN KEY (`idsnpeffeffect`) REFERENCES `snpeffeffect` (`idsnpeffeffect`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `snpeff_snv` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=43134239 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snpeffeffect`
--

DROP TABLE IF EXISTS `snpeffeffect`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snpeffeffect` (
  `idsnpeffeffect` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  `rank` int(11) unsigned DEFAULT NULL,
  PRIMARY KEY (`idsnpeffeffect`)
) ENGINE=InnoDB AUTO_INCREMENT=39 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snv`
--

DROP TABLE IF EXISTS `snv`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snv` (
  `idsnv` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `rs` varchar(15) NOT NULL DEFAULT '',
  `allele` varchar(255) NOT NULL,
  `class` enum('unknown','snp','indel','deletion','cnv') NOT NULL,
  `func` set('unknown','syn','missense','nonsense','stoploss','splice','nearsplice','frameshift','indel','5utr','3utr','noncoding','mirna','lincrna','intronic','intergenic','regulation') NOT NULL,
  `transcript` blob,
  `freq` smallint(6) NOT NULL DEFAULT '0',
  `clinical` tinyint(1) DEFAULT '0',
  `avhet` float NOT NULL DEFAULT '0',
  `valid` set('unknown','by-cluster','by-frequency','by-submitter','by-2hit-2allele','by-hapmap','by-1000genomes') NOT NULL,
  `dp` int(11) NOT NULL DEFAULT '0',
  `af` float NOT NULL DEFAULT '0',
  `refallele` varchar(255) NOT NULL DEFAULT '',
  `lincrna` blob,
  `refseq` blob,
  `length` int(11) unsigned DEFAULT NULL,
  PRIMARY KEY (`idsnv`),
  UNIQUE KEY `uniquesnv` (`chrom`,`start`,`end`,`class`,`refallele`,`allele`),
  KEY `chromstart` (`chrom`,`start`),
  KEY `chromend` (`chrom`,`end`),
  KEY `rs` (`rs`),
  KEY `func` (`func`),
  KEY `class` (`class`),
  KEY `af` (`af`),
  KEY `avhet` (`avhet`)
) ENGINE=InnoDB AUTO_INCREMENT=10990357 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snv2diseasegroup`
--

DROP TABLE IF EXISTS `snv2diseasegroup`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snv2diseasegroup` (
  `idsnv2diseasegroup` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `fidsnv` int(11) unsigned NOT NULL,
  `fiddiseasegroup` int(11) unsigned NOT NULL,
  `fsample` int(11) NOT NULL DEFAULT '0',
  `fsampleall` int(11) NOT NULL DEFAULT '0',
  `fpedigree` int(11) NOT NULL DEFAULT '0',
  `fpedigreeall` int(11) NOT NULL DEFAULT '0',
  `samplecontrols` int(11) NOT NULL DEFAULT '0',
  `pedigreecontrols` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`idsnv2diseasegroup`),
  UNIQUE KEY `diseasegroupsnv` (`fiddiseasegroup`,`fidsnv`),
  KEY `samplecontrolsgroup` (`fiddiseasegroup`,`samplecontrols`,`fidsnv`),
  KEY `pedigreecontrolsgroup` (`fiddiseasegroup`,`pedigreecontrols`,`fidsnv`),
  KEY `s2didsnvgroup` (`fidsnv`),
  KEY `s2dcontrols` (`samplecontrols`),
  CONSTRAINT `s2diddiseasegroup` FOREIGN KEY (`fiddiseasegroup`) REFERENCES `exomehg19`.`diseasegroup` (`iddiseasegroup`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `s2didsnvgroup` FOREIGN KEY (`fidsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=4415847 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snvgene`
--

DROP TABLE IF EXISTS `snvgene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snvgene` (
  `idsnvgene` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idgene` int(11) unsigned NOT NULL,
  `idsnv` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idsnvgene`),
  UNIQUE KEY `snvgene` (`idgene`,`idsnv`),
  KEY `gene` (`idgene`),
  KEY `snv` (`idsnv`),
  KEY `snvgene2` (`idsnv`,`idgene`),
  CONSTRAINT `gene` FOREIGN KEY (`idgene`) REFERENCES `gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `snv` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=12982128 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snvsample`
--

DROP TABLE IF EXISTS `snvsample`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snvsample` (
  `idsnvsample` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsnv` int(11) unsigned NOT NULL,
  `idsample` int(11) unsigned NOT NULL,
  `alleles` tinyint(3) unsigned NOT NULL,
  `percentvar` smallint(5) unsigned NOT NULL,
  `percentfor` float NOT NULL,
  `medianqual` smallint(6) NOT NULL,
  `refqual` smallint(5) unsigned NOT NULL,
  `snvqual` smallint(5) unsigned NOT NULL,
  `mapqual` smallint(5) unsigned NOT NULL,
  `coverage` smallint(5) unsigned NOT NULL,
  `filter` set('PASS','VARQ','VARd','VARD','VARa','VARG','VARg','VARP','VARM','VARS','Q20','Q15','Q10','Q5','Q3','TooManyCNVs','VQSRTrancheINDEL99.00to99.90','VQSRTrancheINDEL99.90to100.00','VQSRTrancheSNP99.00to99.90','VQSRTrancheSNP99.90to100.00','QDfilter','MQfilter','FSfilter','MQRankSumfilter','ReadPosRankSumfilter','HaplotypeScorefilter','InbreedingCoefffilter') DEFAULT NULL,
  `gtqual` smallint(5) unsigned NOT NULL,
  `caller` set('samtools','pindel','exomedepth','gatk','cnvnator','lumpy-sv') DEFAULT NULL,
  PRIMARY KEY (`idsnvsample`),
  UNIQUE KEY `idsnvidsample` (`idsnv`,`idsample`),
  UNIQUE KEY `idsampleidsnv` (`idsample`,`idsnv`),
  KEY `idsample` (`idsample`),
  KEY `idsnv` (`idsnv`),
  KEY `alleles` (`alleles`),
  KEY `test2` (`filter`,`snvqual`,`mapqual`,`caller`,`alleles`,`percentvar`),
  KEY `caller` (`caller`,`percentvar`,`idsnv`,`idsample`),
  CONSTRAINT `snvsample_idsample` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `snvsample_idsnv` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=17777404 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snvsample_st`
--

DROP TABLE IF EXISTS `snvsample_st`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snvsample_st` (
  `idsnvsample` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsnv` int(11) unsigned NOT NULL,
  `idsample` int(11) unsigned NOT NULL,
  `alleles` tinyint(3) unsigned NOT NULL,
  `percentvar` smallint(5) unsigned NOT NULL,
  `percentfor` float NOT NULL,
  `medianqual` smallint(6) NOT NULL,
  `refqual` smallint(5) unsigned NOT NULL,
  `snvqual` smallint(5) unsigned NOT NULL,
  `mapqual` smallint(5) unsigned NOT NULL,
  `coverage` smallint(5) unsigned NOT NULL,
  `filter` set('PASS','VARQ','VARd','VARD','VARa','VARG','VARg','VARP','VARM','VARS','Q20','Q15','Q10','Q5','Q3','TooManyCNVs') DEFAULT NULL,
  `gtqual` smallint(5) unsigned NOT NULL,
  `caller` set('samtools','pindel','exomedepth','gatk','cnvnator','lumpy-sv') DEFAULT NULL,
  PRIMARY KEY (`idsnvsample`),
  UNIQUE KEY `idsnvidsample_st` (`idsnv`,`idsample`),
  UNIQUE KEY `idsampleidsnv_st` (`idsample`,`idsnv`),
  CONSTRAINT `idsample_st` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `idsnv_st` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=10438095 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snvsample_vs`
--

DROP TABLE IF EXISTS `snvsample_vs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snvsample_vs` (
  `idsnvsample` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsnv` int(11) unsigned NOT NULL,
  `idsample` int(11) unsigned NOT NULL,
  `alleles` tinyint(3) unsigned NOT NULL,
  `percentvar` smallint(5) unsigned NOT NULL,
  `percentfor` float NOT NULL,
  `medianqual` smallint(6) NOT NULL,
  `refqual` smallint(5) unsigned NOT NULL,
  `snvqual` smallint(5) unsigned NOT NULL,
  `mapqual` smallint(5) unsigned NOT NULL,
  `coverage` smallint(5) unsigned NOT NULL,
  `filter` set('PASS','VARQ','VARd','VARD','VARa','VARG','VARg','VARP','VARM','VARS','Q20','Q15','Q10','Q5','Q3','TooManyCNVs') DEFAULT NULL,
  `gtqual` smallint(5) unsigned NOT NULL,
  `caller` set('samtools','pindel','exomedepth','gatk','cnvnator','lumpy-sv') DEFAULT NULL,
  PRIMARY KEY (`idsnvsample`),
  UNIQUE KEY `idsnvidsample_vs` (`idsnv`,`idsample`),
  UNIQUE KEY `idsampleidsnv_vs` (`idsample`,`idsnv`),
  CONSTRAINT `idsample_vs` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `idsnv_vs` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=14088590 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snvsampledenovo`
--

DROP TABLE IF EXISTS `snvsampledenovo`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snvsampledenovo` (
  `idsnvsample` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `mothergtqual` smallint(5) unsigned NOT NULL,
  `mothercoverage` smallint(5) unsigned NOT NULL,
  `motherpercentvar` smallint(5) unsigned NOT NULL,
  `fathergtqual` smallint(5) unsigned NOT NULL,
  `fathercoverage` smallint(5) unsigned NOT NULL,
  `fatherpercentvar` smallint(5) unsigned NOT NULL,
  `posteriorprobability` int(10) unsigned DEFAULT NULL,
  `denovoprobability` float unsigned DEFAULT '0',
  PRIMARY KEY (`idsnvsample`),
  CONSTRAINT `idsnvsample` FOREIGN KEY (`idsnvsample`) REFERENCES `snvsample` (`idsnvsample`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `translocation`
--

DROP TABLE IF EXISTS `translocation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `translocation` (
  `idtranslocation` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsample` int(11) unsigned NOT NULL,
  `chrom1` varchar(45) NOT NULL,
  `pos1` int(11) unsigned NOT NULL,
  `chrom2` varchar(45) NOT NULL,
  `pos2` int(11) unsigned NOT NULL,
  `idgene1` int(11) unsigned NOT NULL,
  `isingene1` tinyint(1) unsigned NOT NULL,
  `idgene2` int(11) unsigned NOT NULL,
  `isingene2` tinyint(1) unsigned NOT NULL,
  `varianttype` enum('balanced','unbalanced','unclassified','inter_insertions') DEFAULT NULL,
  `num_discordant` int(11) unsigned DEFAULT NULL,
  `sc1` int(11) unsigned DEFAULT NULL,
  `sc2` int(11) unsigned DEFAULT NULL,
  PRIMARY KEY (`idtranslocation`),
  UNIQUE KEY `uniquetranslocation` (`idsample`,`chrom1`,`pos1`,`chrom2`,`pos2`),
  KEY `chromstart1` (`chrom1`,`pos1`),
  KEY `chromstart2` (`chrom2`,`pos2`),
  KEY `translocationidgene1` (`idgene1`),
  KEY `translocationidgene2` (`idgene2`),
  CONSTRAINT `translocationidgene1` FOREIGN KEY (`idgene1`) REFERENCES `gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `translocationidgene2` FOREIGN KEY (`idgene2`) REFERENCES `gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `translocationidsample` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `variantstat`
--

DROP TABLE IF EXISTS `variantstat`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variantstat` (
  `idsample` int(11) unsigned NOT NULL,
  `snv` int(11) unsigned DEFAULT NULL,
  `indel` int(11) unsigned DEFAULT NULL,
  `pindel` int(11) unsigned DEFAULT NULL,
  `exomedepth` int(11) unsigned DEFAULT NULL,
  `snvgtqual` float DEFAULT NULL,
  `snvdepth` float DEFAULT NULL,
  `inserttime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `updatetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`idsample`),
  CONSTRAINT `idsamplevariantstat` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `variantstat_st`
--

DROP TABLE IF EXISTS `variantstat_st`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variantstat_st` (
  `idsample` int(11) unsigned NOT NULL,
  `snv` int(11) unsigned DEFAULT NULL,
  `indel` int(11) unsigned DEFAULT NULL,
  `pindel` int(11) unsigned DEFAULT NULL,
  `exomedepth` int(11) unsigned DEFAULT NULL,
  `snvgtqual` float DEFAULT NULL,
  `snvdepth` float DEFAULT NULL,
  `inserttime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `updatetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`idsample`),
  CONSTRAINT `idsamplevariantstat_st` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `variantstat_vs`
--

DROP TABLE IF EXISTS `variantstat_vs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `variantstat_vs` (
  `idsample` int(11) unsigned NOT NULL,
  `snv` int(11) unsigned DEFAULT NULL,
  `indel` int(11) unsigned DEFAULT NULL,
  `pindel` int(11) unsigned DEFAULT NULL,
  `exomedepth` int(11) unsigned DEFAULT NULL,
  `snvgtqual` float DEFAULT NULL,
  `snvdepth` float DEFAULT NULL,
  `inserttime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `updatetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`idsample`),
  CONSTRAINT `idsamplevariantstat_vs` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2019-11-25 13:57:15
