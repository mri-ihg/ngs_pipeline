-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: genomegatk
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
) ENGINE=InnoDB AUTO_INCREMENT=5289226 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=44827 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=17 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=3563770 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=882853 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `oldsv`
--

DROP TABLE IF EXISTS `oldsv`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `oldsv` (
  `idsv` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `svtype` enum('DEL','DUP','INS','INV','CNV','mCNV') NOT NULL,
  `svlen` int(11) NOT NULL,
  `af1KG` float unsigned DEFAULT NULL COMMENT 'summed allele frequency of 1000 genome variants overlapping >60%',
  `num1KG` int(10) unsigned DEFAULT NULL COMMENT 'number of 1000 genome variants overlapping >60%',
  `type1KG` set('CNV','DEL','DEL_ALU','DEL_HERV','DEL_LINE1','DEL_SVA','DUP','INS','INV') DEFAULT NULL COMMENT 'types of 1000 genome variants overlapping >60%',
  `freq` smallint(6) NOT NULL DEFAULT '0',
  `giaboverlap` float unsigned DEFAULT '0',
  `lowcomploverlap` float unsigned DEFAULT '0',
  `dgvoverlap` int(10) unsigned DEFAULT '0',
  `gapoverlap` float unsigned DEFAULT '0',
  `gensupdupsoverlap` int(10) unsigned DEFAULT '0',
  `overlaps` set('Promoter','Coding_Region') DEFAULT NULL,
  `afdbsnp` float unsigned DEFAULT NULL COMMENT 'minimum minor allele frequency of dbSNP variants overlapping >60%; 999 if overlap but no AF defined',
  PRIMARY KEY (`idsv`),
  UNIQUE KEY `uniquesv` (`chrom`,`start`,`end`,`svtype`),
  KEY `chromstart` (`chrom`,`start`),
  KEY `chromend` (`chrom`,`end`)
) ENGINE=InnoDB AUTO_INCREMENT=2897741 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `oldsvgene`
--

DROP TABLE IF EXISTS `oldsvgene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `oldsvgene` (
  `idsvgene` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idgene` int(11) unsigned NOT NULL,
  `idsv` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idsvgene`),
  UNIQUE KEY `svgene` (`idgene`,`idsv`),
  KEY `gene` (`idgene`),
  KEY `sv` (`idsv`),
  KEY `svgene2` (`idsv`,`idgene`),
  CONSTRAINT `svgene.gene` FOREIGN KEY (`idgene`) REFERENCES `gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `svgene.sv` FOREIGN KEY (`idsv`) REFERENCES `oldsv` (`idsv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=70453682 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `oldsvsample`
--

DROP TABLE IF EXISTS `oldsvsample`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `oldsvsample` (
  `idsvsample` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsv` int(11) unsigned DEFAULT NULL,
  `idsample` int(11) unsigned NOT NULL,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `svtype` enum('DEL','DUP','INS','INV','CNV','mCNV') NOT NULL,
  `svlen` int(11) NOT NULL,
  `caller` set('samtools','pindel','exomedepth','gatk','cnvnator','lumpy-sv','delly','breakdancer','manta') DEFAULT NULL,
  `edexonnumber` int(11) unsigned DEFAULT NULL COMMENT 'exomedepth: number of exons affected by cnv',
  `edratio` float DEFAULT NULL COMMENT 'exomedepth: ratio of cnv: <1=deletion >1=duplication',
  `cndosage` float unsigned DEFAULT NULL COMMENT 'cnvnator: dosage',
  `cnscore1` double DEFAULT NULL COMMENT 'cnvnator: score1',
  `cnscore2` double DEFAULT NULL COMMENT 'cnvnator: score2',
  `cnscore3` double DEFAULT NULL COMMENT 'cnvnator: score3',
  `cnscore4` double DEFAULT NULL,
  `pialleles` int(10) unsigned DEFAULT NULL COMMENT 'pindel: number of variant alleles',
  `pidp` int(10) unsigned DEFAULT NULL COMMENT 'pindel: number of reads at SV position',
  `pipercentvar` float unsigned DEFAULT NULL COMMENT 'pindel: percent of reads showing variant',
  `bddp` int(10) unsigned DEFAULT NULL COMMENT 'breakdancer: number of reads at SV position',
  `bdor1` varchar(100) DEFAULT NULL COMMENT 'breakdancer: read orientation at breakpoint1',
  `bdor2` varchar(100) DEFAULT NULL COMMENT 'breakdancer: read orientation at breakpoint2',
  `lppe` int(10) unsigned DEFAULT NULL COMMENT 'lumpy-sv: paired-end reads showing evidence',
  `lpsr` int(10) unsigned DEFAULT NULL COMMENT 'lumpy-sv: split reads showing evidence',
  `mtalleles` smallint(5) unsigned DEFAULT NULL COMMENT 'manta: alleles (1-->het,2-->hom)',
  `mtgq` smallint(5) unsigned DEFAULT NULL COMMENT 'manta: genotype quality',
  `cnunique` float DEFAULT NULL,
  PRIMARY KEY (`idsvsample`),
  KEY `idsample` (`idsample`),
  KEY `chromstart` (`idsample`,`chrom`,`start`),
  KEY `chromend` (`idsample`,`chrom`,`end`),
  KEY `uniquesvsample` (`idsample`,`chrom`,`start`,`end`,`svtype`),
  KEY `sv_idsv` (`idsv`),
  KEY `svssvtypechromstartend` (`svtype`,`chrom`,`start`,`end`),
  KEY `svssvtypechromstart` (`svtype`,`chrom`,`start`),
  KEY `svssvtypechromend` (`svtype`,`chrom`,`end`),
  CONSTRAINT `sv_idsample` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `sv_idsv` FOREIGN KEY (`idsv`) REFERENCES `oldsv` (`idsv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=24005377 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=136404467 DEFAULT CHARSET=latin1;
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
  `giab` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`idsnv`),
  UNIQUE KEY `uniquesnv` (`chrom`,`start`,`end`,`class`,`refallele`,`allele`),
  KEY `chromstart` (`chrom`,`start`),
  KEY `chromend` (`chrom`,`end`),
  KEY `rs` (`rs`),
  KEY `func` (`func`),
  KEY `class` (`class`),
  KEY `af` (`af`),
  KEY `avhet` (`avhet`)
) ENGINE=InnoDB AUTO_INCREMENT=96022759 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=136518331 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=187586324 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `snvgenebackup`
--

DROP TABLE IF EXISTS `snvgenebackup`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `snvgenebackup` (
  `idsnvgene` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idgene` int(11) unsigned NOT NULL,
  `idsnv` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idsnvgene`),
  UNIQUE KEY `snvgene` (`idgene`,`idsnv`),
  KEY `gene` (`idgene`),
  KEY `snv` (`idsnv`),
  KEY `snvgene2` (`idsnv`,`idgene`)
) ENGINE=InnoDB AUTO_INCREMENT=118131892 DEFAULT CHARSET=latin1;
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
  KEY `filter` (`idsnv`,`filter`,`mapqual`,`gtqual`),
  CONSTRAINT `idsample` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `idsnv` FOREIGN KEY (`idsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=3408925517 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `sv`
--

DROP TABLE IF EXISTS `sv`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `sv` (
  `idsv` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `svtype` enum('DEL','DUP','INS','INV','CNV','mCNV') NOT NULL,
  `svlen` int(11) NOT NULL,
  `af1KG` float unsigned DEFAULT NULL COMMENT 'summed allele frequency of 1000 genome variants overlapping >60%',
  `num1KG` int(10) unsigned DEFAULT NULL COMMENT 'number of 1000 genome variants overlapping >60%',
  `type1KG` set('CNV','DEL','DEL_ALU','DEL_HERV','DEL_LINE1','DEL_SVA','DUP','INS','INV') DEFAULT NULL COMMENT 'types of 1000 genome variants overlapping >60%',
  `freq` smallint(6) NOT NULL DEFAULT '0',
  `giaboverlap` float unsigned DEFAULT '0',
  `lowcomploverlap` float unsigned DEFAULT '0',
  `dgvoverlap` int(10) unsigned DEFAULT '0',
  `gapoverlap` float unsigned DEFAULT '0',
  `gensupdupsoverlap` int(10) unsigned DEFAULT '0',
  `overlaps` set('Promoter','Coding_Region') DEFAULT NULL,
  `afdbsnp` float unsigned DEFAULT NULL COMMENT 'minimum minor allele frequency of dbSNP variants overlapping >60%; 999 if overlap but no AF defined',
  PRIMARY KEY (`idsv`),
  UNIQUE KEY `uniquesv` (`chrom`,`start`,`end`,`svtype`,`svlen`),
  KEY `chromstart` (`chrom`,`start`),
  KEY `chromend` (`chrom`,`end`)
) ENGINE=InnoDB AUTO_INCREMENT=4426462 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `svgene`
--

DROP TABLE IF EXISTS `svgene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `svgene` (
  `idsvgene` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idgene` int(11) unsigned NOT NULL,
  `idsv` int(11) unsigned NOT NULL,
  PRIMARY KEY (`idsvgene`),
  UNIQUE KEY `svgene` (`idgene`,`idsv`),
  KEY `gene` (`idgene`),
  KEY `sv` (`idsv`),
  KEY `svgene2` (`idsv`,`idgene`)
) ENGINE=InnoDB AUTO_INCREMENT=53432111 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `svsample`
--

DROP TABLE IF EXISTS `svsample`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `svsample` (
  `idsvsample` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `idsv` int(11) unsigned DEFAULT NULL,
  `idsample` int(11) unsigned NOT NULL,
  `chrom` varchar(45) NOT NULL,
  `start` int(11) unsigned NOT NULL,
  `end` int(11) unsigned NOT NULL,
  `svtype` enum('DEL','DUP','INS','INV','CNV','mCNV') NOT NULL,
  `svlen` int(11) NOT NULL,
  `caller` set('samtools','pindel','exomedepth','gatk','cnvnator','lumpy-sv','delly','breakdancer','manta','whamg') DEFAULT NULL,
  `edexonnumber` int(11) unsigned DEFAULT NULL COMMENT 'exomedepth: number of exons affected by cnv',
  `edratio` float DEFAULT NULL COMMENT 'exomedepth: ratio of cnv: <1=deletion >1=duplication',
  `cndosage` float unsigned DEFAULT NULL COMMENT 'cnvnator: dosage',
  `cnscore1` double DEFAULT NULL COMMENT 'cnvnator: score1',
  `cnscore2` double DEFAULT NULL COMMENT 'cnvnator: score2',
  `cnscore3` double DEFAULT NULL COMMENT 'cnvnator: score3',
  `cnscore4` double DEFAULT NULL,
  `pialleles` int(10) unsigned DEFAULT NULL COMMENT 'pindel: number of variant alleles',
  `pidp` int(10) unsigned DEFAULT NULL COMMENT 'pindel: number of reads at SV position',
  `pipercentvar` float unsigned DEFAULT NULL COMMENT 'pindel: percent of reads showing variant',
  `bddp` int(10) unsigned DEFAULT NULL COMMENT 'breakdancer: number of reads at SV position',
  `bdor1` varchar(100) DEFAULT NULL COMMENT 'breakdancer: read orientation at breakpoint1',
  `bdor2` varchar(100) DEFAULT NULL COMMENT 'breakdancer: read orientation at breakpoint2',
  `lppe` int(10) unsigned DEFAULT NULL COMMENT 'lumpy-sv: paired-end reads showing evidence',
  `lpsr` int(10) unsigned DEFAULT NULL COMMENT 'lumpy-sv: split reads showing evidence',
  `mtalleles` smallint(5) unsigned DEFAULT NULL COMMENT 'manta: alleles (1-->het,2-->hom)',
  `mtgq` smallint(5) unsigned DEFAULT NULL COMMENT 'manta: genotype quality',
  `cnunique` float DEFAULT NULL,
  PRIMARY KEY (`idsvsample`),
  KEY `idsample` (`idsample`),
  KEY `chromstart` (`idsample`,`chrom`,`start`),
  KEY `chromend` (`idsample`,`chrom`,`end`),
  KEY `uniquesvsample` (`idsample`,`chrom`,`start`,`end`,`svtype`),
  KEY `sv_idsv` (`idsv`),
  KEY `svssvtypechromstartend` (`svtype`,`chrom`,`start`,`end`),
  KEY `svssvtypechromstart` (`svtype`,`chrom`,`start`),
  KEY `svssvtypechromend` (`svtype`,`chrom`,`end`)
) ENGINE=InnoDB AUTO_INCREMENT=13736040 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=3247 DEFAULT CHARSET=latin1;
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
  `snvqual` float DEFAULT NULL,
  `snvgtqual` float DEFAULT NULL,
  `snvdepth` float DEFAULT NULL,
  `inserttime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `updatetime` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  PRIMARY KEY (`idsample`),
  CONSTRAINT `idsamplevariantstat` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
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

-- Dump completed on 2019-11-25 13:57:14
