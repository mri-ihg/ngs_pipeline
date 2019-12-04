-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: rnamm10
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
  PRIMARY KEY (`idgenebased`),
  UNIQUE KEY `idgeneidsample` (`idgene`,`idsample`),
  UNIQUE KEY `idsampleidgene` (`idsample`,`idgene`),
  KEY `idsample_index` (`idsample`) USING BTREE,
  KEY `idgene_index` (`idgene`) USING BTREE,
  CONSTRAINT `idgene` FOREIGN KEY (`idgene`) REFERENCES `exomemm10`.`gene` (`idgene`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `idsample` FOREIGN KEY (`idsample`) REFERENCES `exomehg19`.`sample` (`idsample`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB AUTO_INCREMENT=12806535 DEFAULT CHARSET=utf8;
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
