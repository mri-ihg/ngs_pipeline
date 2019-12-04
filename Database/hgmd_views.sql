-- MySQL dump 10.15  Distrib 10.0.35-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: hgmd_views
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
-- Temporary table structure for view `data_for_ngs`
--

DROP TABLE IF EXISTS `data_for_ngs`;
/*!50001 DROP VIEW IF EXISTS `data_for_ngs`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `data_for_ngs` (
  `hgmd_accession` tinyint NOT NULL,
  `gene_symbol` tinyint NOT NULL,
  `gene_strand` tinyint NOT NULL,
  `chromosome` tinyint NOT NULL,
  `hg19_start` tinyint NOT NULL,
  `hg19_end` tinyint NOT NULL,
  `hgmd_ref` tinyint NOT NULL,
  `hgmd_alt` tinyint NOT NULL,
  `refseq` tinyint NOT NULL,
  `hgvs` tinyint NOT NULL,
  `variant_class` tinyint NOT NULL,
  `dbSNP` tinyint NOT NULL,
  `primary_pubmed` tinyint NOT NULL,
  `additional_pubmed` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `gene_to_concept`
--

DROP TABLE IF EXISTS `gene_to_concept`;
/*!50001 DROP VIEW IF EXISTS `gene_to_concept`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `gene_to_concept` (
  `gene_sym` tinyint NOT NULL,
  `phenotype` tinyint NOT NULL,
  `cui` tinyint NOT NULL,
  `str` tinyint NOT NULL,
  `tty` tinyint NOT NULL,
  `sab` tinyint NOT NULL,
  `code` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `isoform_list`
--

DROP TABLE IF EXISTS `isoform_list`;
/*!50001 DROP VIEW IF EXISTS `isoform_list`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `isoform_list` (
  `hgmd_gene_id` tinyint NOT NULL,
  `gene_symbol` tinyint NOT NULL,
  `gene_description` tinyint NOT NULL,
  `entrezID` tinyint NOT NULL,
  `hgmd_accession` tinyint NOT NULL,
  `variant_class` tinyint NOT NULL,
  `mutation_description` tinyint NOT NULL,
  `refseq` tinyint NOT NULL,
  `hgvs` tinyint NOT NULL,
  `chromosome` tinyint NOT NULL,
  `hg19_start` tinyint NOT NULL,
  `hg19_end` tinyint NOT NULL,
  `primary_pubmed` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `mut_to_concept`
--

DROP TABLE IF EXISTS `mut_to_concept`;
/*!50001 DROP VIEW IF EXISTS `mut_to_concept`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `mut_to_concept` (
  `acc_num` tinyint NOT NULL,
  `gene_sym` tinyint NOT NULL,
  `phenotype` tinyint NOT NULL,
  `rela` tinyint NOT NULL,
  `cui` tinyint NOT NULL,
  `str` tinyint NOT NULL,
  `ispref` tinyint NOT NULL,
  `sab` tinyint NOT NULL,
  `code` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `mut_to_no_cui`
--

DROP TABLE IF EXISTS `mut_to_no_cui`;
/*!50001 DROP VIEW IF EXISTS `mut_to_no_cui`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `mut_to_no_cui` (
  `acc_num` tinyint NOT NULL,
  `gene_sym` tinyint NOT NULL,
  `phenotype` tinyint NOT NULL,
  `rela` tinyint NOT NULL,
  `cui` tinyint NOT NULL,
  `str` tinyint NOT NULL,
  `ispref` tinyint NOT NULL,
  `sab` tinyint NOT NULL,
  `code` tinyint NOT NULL,
  `lat` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `with_additional_references`
--

DROP TABLE IF EXISTS `with_additional_references`;
/*!50001 DROP VIEW IF EXISTS `with_additional_references`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `with_additional_references` (
  `hgmd_accession` tinyint NOT NULL,
  `primary_phenotype` tinyint NOT NULL,
  `additional_phenotype` tinyint NOT NULL,
  `gene_symbol` tinyint NOT NULL,
  `mutation_description` tinyint NOT NULL,
  `refseq` tinyint NOT NULL,
  `hgvs` tinyint NOT NULL,
  `variant_class` tinyint NOT NULL,
  `primary_pubmed` tinyint NOT NULL,
  `additional_pubmed` tinyint NOT NULL,
  `additional_reference_types` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `with_dbsnp_rs`
--

DROP TABLE IF EXISTS `with_dbsnp_rs`;
/*!50001 DROP VIEW IF EXISTS `with_dbsnp_rs`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `with_dbsnp_rs` (
  `hgmd_accession` tinyint NOT NULL,
  `gene_symbol` tinyint NOT NULL,
  `mutation_description` tinyint NOT NULL,
  `refseq` tinyint NOT NULL,
  `hgvs` tinyint NOT NULL,
  `variant_class` tinyint NOT NULL,
  `dbSNP` tinyint NOT NULL,
  `1000G_frequency` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Temporary table structure for view `with_edit_history`
--

DROP TABLE IF EXISTS `with_edit_history`;
/*!50001 DROP VIEW IF EXISTS `with_edit_history`*/;
SET @saved_cs_client     = @@character_set_client;
SET character_set_client = utf8;
/*!50001 CREATE TABLE `with_edit_history` (
  `gene_symbol` tinyint NOT NULL,
  `hgmd_accession` tinyint NOT NULL,
  `mutation_description` tinyint NOT NULL,
  `hgvs` tinyint NOT NULL,
  `variant_class` tinyint NOT NULL,
  `column_updated` tinyint NOT NULL,
  `before_update` tinyint NOT NULL,
  `after_update` tinyint NOT NULL,
  `date_updated` tinyint NOT NULL
) ENGINE=MyISAM */;
SET character_set_client = @saved_cs_client;

--
-- Final view structure for view `data_for_ngs`
--

/*!50001 DROP TABLE IF EXISTS `data_for_ngs`*/;
/*!50001 DROP VIEW IF EXISTS `data_for_ngs`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `data_for_ngs` AS select `hgmd_pro`.`hg19_coords`.`acc_num` AS `hgmd_accession`,`hgmd_pro`.`allmut`.`gene` AS `gene_symbol`,`hgmd_pro`.`hg19_coords`.`strand` AS `gene_strand`,`hgmd_pro`.`hg19_coords`.`chromosome` AS `chromosome`,`hgmd_pro`.`hg19_coords`.`coordSTART` AS `hg19_start`,`hgmd_pro`.`hg19_coords`.`coordEND` AS `hg19_end`,ucase(`hgmd_pro`.`mutnomen`.`wildBASE`) AS `hgmd_ref`,ucase(`hgmd_pro`.`mutnomen`.`mutBASE`) AS `hgmd_alt`,concat(`hgmd_pro`.`mutnomen`.`refCORE`,'.',`hgmd_pro`.`mutnomen`.`refVER`) AS `refseq`,concat('c.',`hgmd_pro`.`mutnomen`.`hgvs`) AS `hgvs`,`hgmd_pro`.`allmut`.`tag` AS `variant_class`,`hgmd_pro`.`allmut`.`dbsnp` AS `dbSNP`,`hgmd_pro`.`allmut`.`pmid` AS `primary_pubmed`,group_concat(`hgmd_pro`.`extrarefs`.`pmid` separator ', ') AS `additional_pubmed` from (((`hgmd_pro`.`hg19_coords` join `hgmd_pro`.`mutnomen` on((`hgmd_pro`.`hg19_coords`.`acc_num` = `hgmd_pro`.`mutnomen`.`acc_num`))) join `hgmd_pro`.`allmut` on((`hgmd_pro`.`hg19_coords`.`acc_num` = `hgmd_pro`.`allmut`.`acc_num`))) left join `hgmd_pro`.`extrarefs` on((`hgmd_pro`.`hg19_coords`.`acc_num` = `hgmd_pro`.`extrarefs`.`acc_num`))) group by `hgmd_pro`.`hg19_coords`.`acc_num` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `gene_to_concept`
--

/*!50001 DROP TABLE IF EXISTS `gene_to_concept`*/;
/*!50001 DROP VIEW IF EXISTS `gene_to_concept`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `gene_to_concept` AS select `hgmd_phenbase`.`hgmd_mutation`.`gene_sym` AS `gene_sym`,`hgmd_phenbase`.`hgmd_phenotype`.`phenotype` AS `phenotype`,`hgmd_phenbase`.`phenotype_concept`.`cui` AS `cui`,`hgmd_phenbase`.`concept`.`str` AS `str`,`hgmd_phenbase`.`concept`.`tty` AS `tty`,`hgmd_phenbase`.`concept`.`sab` AS `sab`,`hgmd_phenbase`.`concept`.`code` AS `code` from (((`hgmd_phenbase`.`hgmd_mutation` left join `hgmd_phenbase`.`hgmd_phenotype` on((`hgmd_phenbase`.`hgmd_mutation`.`phen_id` = `hgmd_phenbase`.`hgmd_phenotype`.`phen_id`))) left join `hgmd_phenbase`.`phenotype_concept` on((`hgmd_phenbase`.`phenotype_concept`.`phen_id` = `hgmd_phenbase`.`hgmd_mutation`.`phen_id`))) left join `hgmd_phenbase`.`concept` on((`hgmd_phenbase`.`phenotype_concept`.`cui` = `hgmd_phenbase`.`concept`.`cui`))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `isoform_list`
--

/*!50001 DROP TABLE IF EXISTS `isoform_list`*/;
/*!50001 DROP VIEW IF EXISTS `isoform_list`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `isoform_list` AS select ifnull(`hgmd_pro`.`markname`.`svar`,`hgmd_pro`.`markname`.`gene_id`) AS `hgmd_gene_id`,`hgmd_pro`.`markname`.`genesym` AS `gene_symbol`,`hgmd_pro`.`markname`.`genename` AS `gene_description`,`hgmd_pro`.`markname`.`entrezID` AS `entrezID`,`hgmd_pro`.`allmut`.`acc_num` AS `hgmd_accession`,`hgmd_pro`.`allmut`.`tag` AS `variant_class`,`hgmd_pro`.`allmut`.`descr` AS `mutation_description`,concat(`hgmd_pro`.`gene2refseq`.`refcore`,'.',`hgmd_pro`.`gene2refseq`.`refversion`) AS `refseq`,concat('c.',`hgmd_pro`.`allmut`.`hgvs`) AS `hgvs`,`hgmd_pro`.`allmut`.`chromosome` AS `chromosome`,`hgmd_pro`.`allmut`.`startCoord` AS `hg19_start`,`hgmd_pro`.`allmut`.`endCoord` AS `hg19_end`,`hgmd_pro`.`allmut`.`pmid` AS `primary_pubmed` from ((`hgmd_pro`.`markname` join `hgmd_pro`.`allmut`) join `hgmd_pro`.`gene2refseq`) where ((`hgmd_pro`.`markname`.`gene_id` in (select `hgmd_pro`.`markname`.`svar` from `hgmd_pro`.`markname` where (`hgmd_pro`.`markname`.`svar` is not null)) or (`hgmd_pro`.`markname`.`svar` is not null)) and (`hgmd_pro`.`markname`.`genesym` = `hgmd_pro`.`allmut`.`gene`) and (`hgmd_pro`.`markname`.`gene_id` = `hgmd_pro`.`gene2refseq`.`hgmdID`)) order by `hgmd_pro`.`markname`.`genesym` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `mut_to_concept`
--

/*!50001 DROP TABLE IF EXISTS `mut_to_concept`*/;
/*!50001 DROP VIEW IF EXISTS `mut_to_concept`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `mut_to_concept` AS select `hgmd_phenbase`.`hgmd_mutation`.`acc_num` AS `acc_num`,`hgmd_phenbase`.`hgmd_mutation`.`gene_sym` AS `gene_sym`,`hgmd_phenbase`.`hgmd_phenotype`.`phenotype` AS `phenotype`,`hgmd_phenbase`.`phenotype_concept`.`rela` AS `rela`,`hgmd_phenbase`.`phenotype_concept`.`cui` AS `cui`,`hgmd_phenbase`.`concept`.`str` AS `str`,`hgmd_phenbase`.`concept`.`ispref` AS `ispref`,`hgmd_phenbase`.`concept`.`sab` AS `sab`,`hgmd_phenbase`.`concept`.`code` AS `code` from (((`hgmd_phenbase`.`hgmd_mutation` join `hgmd_phenbase`.`hgmd_phenotype` on((`hgmd_phenbase`.`hgmd_mutation`.`phen_id` = `hgmd_phenbase`.`hgmd_phenotype`.`phen_id`))) join `hgmd_phenbase`.`phenotype_concept` on((`hgmd_phenbase`.`phenotype_concept`.`phen_id` = `hgmd_phenbase`.`hgmd_mutation`.`phen_id`))) join `hgmd_phenbase`.`concept` on((`hgmd_phenbase`.`phenotype_concept`.`cui` = `hgmd_phenbase`.`concept`.`cui`))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `mut_to_no_cui`
--

/*!50001 DROP TABLE IF EXISTS `mut_to_no_cui`*/;
/*!50001 DROP VIEW IF EXISTS `mut_to_no_cui`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `mut_to_no_cui` AS select `hgmd_phenbase`.`hgmd_mutation`.`acc_num` AS `acc_num`,`hgmd_phenbase`.`hgmd_mutation`.`gene_sym` AS `gene_sym`,`hgmd_phenbase`.`hgmd_phenotype`.`phenotype` AS `phenotype`,`hgmd_phenbase`.`phenotype_concept`.`rela` AS `rela`,`hgmd_phenbase`.`phenotype_concept`.`cui` AS `cui`,`hgmd_phenbase`.`concept`.`str` AS `str`,`hgmd_phenbase`.`concept`.`ispref` AS `ispref`,`hgmd_phenbase`.`concept`.`sab` AS `sab`,`hgmd_phenbase`.`concept`.`code` AS `code`,`hgmd_phenbase`.`concept`.`lat` AS `lat` from (((`hgmd_phenbase`.`hgmd_mutation` left join `hgmd_phenbase`.`hgmd_phenotype` on((`hgmd_phenbase`.`hgmd_mutation`.`phen_id` = `hgmd_phenbase`.`hgmd_phenotype`.`phen_id`))) left join `hgmd_phenbase`.`phenotype_concept` on((`hgmd_phenbase`.`phenotype_concept`.`phen_id` = `hgmd_phenbase`.`hgmd_mutation`.`phen_id`))) left join `hgmd_phenbase`.`concept` on((`hgmd_phenbase`.`phenotype_concept`.`cui` = `hgmd_phenbase`.`concept`.`cui`))) where isnull(`hgmd_phenbase`.`phenotype_concept`.`cui`) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `with_additional_references`
--

/*!50001 DROP TABLE IF EXISTS `with_additional_references`*/;
/*!50001 DROP VIEW IF EXISTS `with_additional_references`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `with_additional_references` AS select `hgmd_pro`.`allmut`.`acc_num` AS `hgmd_accession`,`hgmd_pro`.`allmut`.`disease` AS `primary_phenotype`,group_concat(distinct `hgmd_pro`.`extrarefs`.`disease` separator '; ') AS `additional_phenotype`,`hgmd_pro`.`allmut`.`gene` AS `gene_symbol`,`hgmd_pro`.`allmut`.`descr` AS `mutation_description`,concat(`hgmd_pro`.`mutnomen`.`refCORE`,'.',`hgmd_pro`.`mutnomen`.`refVER`) AS `refseq`,concat('c.',`hgmd_pro`.`mutnomen`.`hgvs`) AS `hgvs`,`hgmd_pro`.`allmut`.`tag` AS `variant_class`,`hgmd_pro`.`allmut`.`pmid` AS `primary_pubmed`,group_concat(distinct `hgmd_pro`.`extrarefs`.`pmid` separator '; ') AS `additional_pubmed`,group_concat(distinct `hgmd_pro`.`extrarefs`.`reftag` separator '; ') AS `additional_reference_types` from ((`hgmd_pro`.`allmut` join `hgmd_pro`.`mutnomen` on((`hgmd_pro`.`allmut`.`acc_num` = `hgmd_pro`.`mutnomen`.`acc_num`))) join `hgmd_pro`.`extrarefs` on((`hgmd_pro`.`allmut`.`acc_num` = `hgmd_pro`.`extrarefs`.`acc_num`))) group by `hgmd_pro`.`allmut`.`acc_num` */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `with_dbsnp_rs`
--

/*!50001 DROP TABLE IF EXISTS `with_dbsnp_rs`*/;
/*!50001 DROP VIEW IF EXISTS `with_dbsnp_rs`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `with_dbsnp_rs` AS select `hgmd_pro`.`allmut`.`acc_num` AS `hgmd_accession`,`hgmd_pro`.`allmut`.`gene` AS `gene_symbol`,`hgmd_pro`.`allmut`.`descr` AS `mutation_description`,concat(`hgmd_pro`.`mutnomen`.`refCORE`,'.',`hgmd_pro`.`mutnomen`.`refVER`) AS `refseq`,concat('c.',`hgmd_pro`.`mutnomen`.`hgvs`) AS `hgvs`,`hgmd_pro`.`allmut`.`tag` AS `variant_class`,`hgmd_pro`.`dbsnp`.`dbsnp_id` AS `dbSNP`,`hgmd_pro`.`dbsnp`.`MAFfreq` AS `1000G_frequency` from ((`hgmd_pro`.`allmut` left join `hgmd_pro`.`mutnomen` on((`hgmd_pro`.`allmut`.`acc_num` = `hgmd_pro`.`mutnomen`.`acc_num`))) join `hgmd_pro`.`dbsnp` on((`hgmd_pro`.`allmut`.`acc_num` = `hgmd_pro`.`dbsnp`.`hgmd_acc`))) */;
/*!50001 SET character_set_client      = @saved_cs_client */;
/*!50001 SET character_set_results     = @saved_cs_results */;
/*!50001 SET collation_connection      = @saved_col_connection */;

--
-- Final view structure for view `with_edit_history`
--

/*!50001 DROP TABLE IF EXISTS `with_edit_history`*/;
/*!50001 DROP VIEW IF EXISTS `with_edit_history`*/;
/*!50001 SET @saved_cs_client          = @@character_set_client */;
/*!50001 SET @saved_cs_results         = @@character_set_results */;
/*!50001 SET @saved_col_connection     = @@collation_connection */;
/*!50001 SET character_set_client      = utf8 */;
/*!50001 SET character_set_results     = utf8 */;
/*!50001 SET collation_connection      = utf8_general_ci */;
/*!50001 CREATE ALGORITHM=UNDEFINED */
/*!50013 DEFINER=`root`@`localhost` SQL SECURITY DEFINER */
/*!50001 VIEW `with_edit_history` AS select `hgmd_pro`.`allmut`.`gene` AS `gene_symbol`,`hgmd_pro`.`allmut`.`acc_num` AS `hgmd_accession`,`hgmd_pro`.`allmut`.`descr` AS `mutation_description`,concat('c.',`hgmd_pro`.`allmut`.`hgvs`) AS `hgvs`,`hgmd_pro`.`allmut`.`tag` AS `variant_class`,`hgmd_pro`.`history`.`colname` AS `column_updated`,`hgmd_pro`.`history`.`beforeUpd` AS `before_update`,`hgmd_pro`.`history`.`afterUpd` AS `after_update`,`hgmd_pro`.`history`.`updDate` AS `date_updated` from (`hgmd_pro`.`allmut` join `hgmd_pro`.`history`) where (`hgmd_pro`.`allmut`.`acc_num` = `hgmd_pro`.`history`.`acc_num`) order by `hgmd_pro`.`allmut`.`gene`,`hgmd_pro`.`allmut`.`acc_num`,`hgmd_pro`.`history`.`updDate` */;
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
