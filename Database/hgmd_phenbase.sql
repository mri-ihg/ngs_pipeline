-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: hgmd_phenbase
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
-- Table structure for table `concept`
--

DROP TABLE IF EXISTS `concept`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `concept` (
  `cui` char(10) COLLATE utf8_unicode_ci DEFAULT NULL,
  `str` text COLLATE utf8_unicode_ci,
  `sab` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `code` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `aui` varchar(15) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `ispref` char(1) COLLATE utf8_unicode_ci DEFAULT NULL,
  `ts` char(1) COLLATE utf8_unicode_ci NOT NULL,
  `tty` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `stt` varchar(3) COLLATE utf8_unicode_ci NOT NULL,
  `lat` char(3) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`aui`),
  KEY `cui` (`cui`),
  KEY `sab` (`sab`),
  KEY `code` (`code`),
  KEY `str` (`str`(30)),
  KEY `tty` (`tty`),
  KEY `ts` (`ts`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `concept_def`
--

DROP TABLE IF EXISTS `concept_def`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `concept_def` (
  `cui` char(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `sab` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `def` text COLLATE utf8_unicode_ci,
  `aui` varchar(9) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `atui` varchar(11) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `satui` varchar(50) COLLATE utf8_unicode_ci DEFAULT NULL,
  `suppress` char(1) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`atui`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hgmd_mutation`
--

DROP TABLE IF EXISTS `hgmd_mutation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hgmd_mutation` (
  `acc_num` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  `gene_sym` varchar(15) COLLATE utf8_unicode_ci DEFAULT NULL,
  `phen_id` int(11) NOT NULL,
  PRIMARY KEY (`acc_num`,`phen_id`),
  KEY `acc_num` (`acc_num`),
  KEY `gene_sym` (`gene_sym`),
  KEY `hgmd_disease` (`phen_id`),
  KEY `acc_num_2` (`acc_num`,`gene_sym`,`phen_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `hgmd_phenotype`
--

DROP TABLE IF EXISTS `hgmd_phenotype`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `hgmd_phenotype` (
  `phen_id` int(11) NOT NULL,
  `phenotype` varchar(125) COLLATE utf8_unicode_ci DEFAULT NULL,
  PRIMARY KEY (`phen_id`),
  KEY `k_hgmd_phenotype_concept` (`phenotype`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `phenotype_concept`
--

DROP TABLE IF EXISTS `phenotype_concept`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `phenotype_concept` (
  `phen_id` int(11) NOT NULL,
  `rela` enum('root','is_a','is_not','sibling_of') COLLATE utf8_unicode_ci NOT NULL,
  `cui` varchar(10) COLLATE utf8_unicode_ci NOT NULL DEFAULT '',
  PRIMARY KEY (`phen_id`,`rela`,`cui`),
  KEY `my_cui_key` (`cui`)
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
