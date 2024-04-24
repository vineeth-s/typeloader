# Version History

## 2.15.0 (2024-04-24)
## Features:
- added support for DRB345 [TL2-286]

## Bugfixes:
- add warning for 5' overhang edgecase [TL2-293]
- unnecessary popup warning for LSL-internal users removed [TL2-289]

## Backend:
- adjusted to changes of hla.dat by IPD [TL2-292] 
- ENA Webin-CLI update to 7.1.1 [TL2-291]
- poetry update [TL2-296]

## 2.14.2 (2023-10-11)
## Bugfixes:
- retrieve ENA accessions via proxy if specified [TL2-278]
- handle timeout of ENA-request during IPD submission [TL2-282]

### Backend:
- update ENA webin CLI to 6.5.0 [TL2-270]
- poetry update [TL2-283]

## 2.14.1 (2023-06-20)
### Bugfixes:
- enforce correct header in SampleView [TL2-271]

## 2.14.0 (2023-06-02)
### Features:
- retrieve ENA accession numbers directly from server, to compensate for a change at ENA (239)
- explicitly set ENA submissions to publish immediately, to compensate for a change at ENA [TL2-242]
- add handling of spatiotemporal data [TL2-250, TL2-247, TL2-248, TL2-251, TL2-249, TL2-252, TL2-243, TL2-259]
- store reference md5 in text file instead of recomputing it => faster startup (#237)

## Docu:
- updated Linux setup instructions (#234) 

### Backend:
- update ENA webin CLI to 6.4.1 [TL2-256]
- poetry update (#238, TL2-258)

## 2.13.2 (2022-08-31)
### Bugfixes for local:
- set KIR2DL5 pretypings to "unknown" to avoid phasing ambiguities (#230)
- add DQA1 and DPA1 to pretypings (#228)

### Backend:
- update ENA webin CLI to 5.1.0 (#231)

## 2.13.1 (2022-06-29)
- Bugfix for local: handle invalid KIR2DL5 pretypings in the local db (#219)

## 2.13.0 (2022-06-01)
### Backend:
- now uses Python 3.10 (#215)

## 2.12.2 (2022-05-24)
- adjusted required loci for KIR alleles to new IPD standards (#211)

### Backend:
- Updated ENA CLI to newest version 4.4.0 (#214)

## 2.12.1 (2022-02-18)
### Bugfixes:
- fix pretypings files for HLA-DMA, -DMB, -DOA, -DOB (#208)
- fix "all IPD-required loci present in pretypings" check

## 2.12.0 (2022-02-18)
### New Features:
- Added support for non-classical HLA class II genes HLA-DMA, -DMB, -DOA, -DOB (#208)

## Improvements & Bugfixes:
- improve timeout for reference downloads (#205)
  - increase timeout to 60 seconds to help users on a slow connection
  - handle timeout gracefully

### Backend:
- Updated ENA CLI to newest version 4.3.0 (#209)
- route all calls to external sites through proxy if specified (#206) 

## 2.11.1 (2021-11-25)
- Updated ENA CLI to newest version 4.2.3 (#201)

## Improvements:
- rename webin-cli.report after submissions to avoid overwriting (#199)
- Extended logging for webin-cli call (#198)

## Bugfixes:
- if playsound doesn't work, just skip it
- check for java installation before 1st webin-cli call

## 2.11.0 (2021-11-04)
### New Features:
- ENA files now include the TypeLoader & reference database version used to create them (#196)

### Improvements:
- RestartAllele Workflow adjusted to the fact that ENA does not allow sequence updates anymore (#85, #176)
- Throw meaningful errors if:
  - the input XML file cannot be handled (#191)
  - the reference db has probably gone stale (#193)

### Backend:
- Updated ENA CLI to newest version 4.2.1 (#192)
 - change: fusion introns can no longer be submitted to ENA as, e.g., 3/4 and are therefore submitted as 3 (but continue to be submitted as 3/4 to IPD) (#194)
 - submission software & version can now be transmitted to ENA via flatfile (#196)

## 2.10.0 (2021-02-16)
### New Features:
- Option to compare 2 TypeLoader generated files of the same type with each other via SampleView => View Files (#184)

### Improvements:
- Reference db updates between TypeLoader-Upload and IPD Submissions are now handled decently instead of throwing a vague BLAST error (#184)

## 2.9.0 (2020-12-10)
### New Features:
- New option to manually reset reference database to specific version (#171)
- New option to show all currently used versions (for reference databases and TypeLoader itself) (#174) 
- The used reference version is now included in IPD files (#172)
- New "Restart Allele" option. This provides a consistent way to handle necessary sequence corrections whenever the originally used sequence turns out to have been incorrect. By using this option, the ENA and IPD submission numbers are kept but otherwise the allele is started and all files generated from scratch with a fresh input sequence file. The allele gets the comment "restarted with fresh sequence". (#85)
- As a consequence, allele- and project files are no longer editable through TypeLoader. If you need to make edits, use the Restart Allele option.
- Comments are now displayed in the "General" tab of the SampleView insted of the "Lab Processing" tab. This should make it more prominent. 
- New projects: if the optional field "title" is left empty, the content of the field "pool" is copied there now. This enables easier identification of projects in the ENA webportal. (#177)

### Bugfixes:
- allele numbering after allele deletions has been fixed (#181)

### Backend:
- Updated ENA CLI to newest version 3.2.2

## 2.8.3.1 (2020-08-31)
- several small bugfixes

## 2.8.3 (2020-10-27)
- ask for confirmation before submissions to ENA's productive server, to avoid accidental submissions (#170)

## 2.8.2 (2020-08-26)
### New Features:
- Support for ENA Webin-CLI 3.1.0 (#154)
- Timeout after ENA server can't be reached for a certain time. Threshold can be set via Settings => Preferences (#161)
- Whitespaces in pretypings file are now removed (#165)

### Backend:
- Disabled "number of CDS bases must be divisable by 3" check for null alleles for HLA-E, to accomodate for a possible IMGT bug that previously made all HLA-E alleles null-alleles (#162)
- Converted to poetry project 

## 2.8.1 (2020-07-27)
### New Features:
- Added support for non-classical HLA genes HLA-F, -G, -H, -K, -J (#155)
- New option: rebuild reference manually (#156)
- Play an unobtrusive sound when long jobs finish (#158) 
- Support for ENA Webin-CLI 3.0.1 (#154)
- Example files now have a "download all" option (#160)

### Bugfixes:
- BLAST output is now caught in the log (#157)
- Fixed Bug in broken-reference-recognition of KIR

## 2.8.0 (2020-06-10)
### New Features:
- Homozygous XML files are now accepted as input (#148)
- Alleles whose closest known allele is not found automatically due to reference database inconsistencies can now be added via a restricted reference with reference alleles chosen by the user. (#149)

### Bugfix:
- handle unexpected replies from ENA during project creation (#151)

## 2.7.1 (2020-03-09)
### Bugfix: 
- handle Webin-CLI replies with non-standard formatting (#147)

## 2.7.0 (2020-02-27)
### New Features:
- KIR2DL5A and 2DL5B are now treated as separate loci, in accordance with IPD's wishes.
- Ignore KIR pretypings for non-KIR IPD files, because IMGT/HLA cannot handle GL-strings. (#121)
- If a lock is encountered during IPD submission, it is now shown who created the lock. This enables the user to better decide whether this is a crash leftover and can safely be removed. (#136)
- If an allele is too different from all officially known alleles to handle correctly, this is now recognized and to the user. (#138) 
- The 'Options' menu now offers the downloading of logfiles. (#105)

### BugFixes:
- TypeLoader now handles changes within the first or last 3 bp of the allele correctly. (#124)
- Cleaner handling of XML input files: (#115)
  - Apply sanity check to both alleles, not just the first.
  - After user has chosen which allele to submit, delete the other one from this allele's files.
  - Thereby, both alleles of one locus of one sample starting with the same first field are now handled correctly.
- Alleles which start (slightly) before the 5' UTR of the reference allele are now handled correctly. (Previously the codons in the IPD difference strings were off.) (#143)
- Whitespaces in input-filenames can now be handled. (#134)
- Some adjustments for communication with the new ENA Webin CLI version. (#135)
- Safety: Closed projects can no longer be used to add or sumbit alleles. (#58)

## 2.6.1. (2020-01-14)
- BugFix: if rensponse from ENA's WebinCLI cannot be parsed, display it directly. (#135)

## 2.6.0 (2020-01-12)
- Added support for HLA-E (#131) and DPA1 & DQA1 (#118)
- Updated ENA Webin CLI to V2.2.0 (#129)
- Adjusted the use of the tags "gene" and "pseudogene" in ENA files to the requirements of the new Webin CLI. (#130)

## V2.5.1 (2019-10-21)
- During IPD file creation, reject GL strings in non-KIR loci (#123)

## V2.5.0 (2019-09-11)
### Bugfixes:
- in IPD files, now the correct pretyping delimiter (comma or plus) is used even if both alleles are novel (#122) 
- correctly catch and display errors from ENA submission, if they happen

## V2.4.1 (2019-08-07)
### Changes:
- the Windows installer is now provided via GitHub release and is no longer part of the repository itself (#120)

## V2.4.0 (2019-07-15)
### New features:
- ENA submission is now handled through the new Webin-CLI (#100)
- added support for DPA1 and DQA1 (#118)

### Changes & adjustments:
- pretypings are now limited to HLA, KIR and MIC, as per IPD request (#119)
- if allele belongs to an unsupported locus, raise a meaningful error (#118)
- New Projects: all fields are now checked for disallowed characters (#116)

### Behind the scenes:
- moved ENA submission functionality from ENASubmissionForm to typeloader_functions (#14) 

## V2.3.2.1 (2019-06-18)
### Bugfix:
- target alleles with a deletion and a mismatch near the end of the CDS can now be handled (#113)

## V2.3.2 (2019-05-29)
### New features:
- users can delete individual alleles or all alleles of a project by rightclicking in the navigation area (on an allele or project, respectively)
- DKMS-users can now retrieve pretypings directly from the database (#16, #84)
- BothAllelesNovelDialog is now scrollable (#112)

### Changes
- IPD files: fields "name" and "lab\_contact" now contain the general company identifier, while the individual submitter is given as "alt\_contact" (#108)
- all pretypings in IPD files are now consistently formatted as GL-strings
- confirmation IPD files now contain the closest known allele in the description section (#106)

### Bugfixes
- TypeLoader errors during reference updates are now caught (previously led to crashes) (#109)

## V2.3.1.1 (2019-05-07)
### Bugfixes
- various bugfixes to make the Windows installer work

## V2.3.1 (2019-05-03)
### new features
- Format of pretypings file changed (#95)
- IPD reference data: TypeLoader now uses only tested-and-ok versions of IPD's releases, to avoid automatically using broken releases (#99)
- added example files for all gene systems and steps (#88)
- sanity-check for fasta file format (#75)

### Bugfixes:
- the key "sample_id_int" in fasta headers is now recognized correctly (#97)
- handle fasta file extensions beside .fa (#87)
- replaced .locked with _locked to hopefully work for external Windows users (#101)

## 2.3.0 (2019-04-18)
### new features:
- Support for MICA and MICB (#52)
- Support for incomplete sequences: if only the UTRs are partial, TypeLoader now throws a warning but the user can decide to upload anyway. Partially known exons or introns are still rejected (and will continue to be). (#63)

### Bugfixes:
- use class II enumeration for all HLA class II genes, not just DPB1
- many small ones

## 2.2.1 (2019-03-11)
### new features:
- Recognize multiple novel alleles in target locus: this is now checked during IPD file creation amd the files are formatted accordingly. If necessary, user input is requested to clarify which of the given alleles is the target allele. (#72)
- Pretyping checks: TypeLoader now checks whether the data given in the pretypings file fits the target alleles, and requests user changes where necessary to ensure consistent, sensible output. (#82)
- Faster startup: the AlleleOverview, which is used rarely but time-intensive to build, is now only generated on demand. (#80)
- Options for Workflow settings are now consistent with options available through IPD's web form (#78)

### Updated user manual:
- Description for pretypings file is much more elaborate now.

## 2.2.0 (2018-12-19)
### new features:
- automatic detection of null alleles (#6)
- policy change with major structural changes: the local_name is now the unique allele identifier, whereas cell_line is used for samples (#64)
- IPD filenames, cell lines and local_names are now generated automatically (#57, #64)
- added allele history tab to AlleleView (#10)
- added dialog to view user manual (#70)
- download an example pretypings file (#55)
- checks and popups to guide new users (#54) 
 
### Bugfixes:
- sequences with incomplete UTR3 are now rejected (previously, they were accepted but produced buggy outcome) (#62)
- workaround for ENA bug (certain sequences were unfairly rejected) (#53)
	
### Updated user manual:
- new: future features
- new: elaborate first_start page to get new users started
- new: elaborate instructions how to feed TypeLoader-generated data into IPD's web form
- adjusted to all changes 

## 2.1.0 (2018-10-23)
### Bugfixes:
- Codon counting on DPB1 is now correct (#1) 
- added missing files and dirs:
 - Windows installer (#42)
 - table .csv files (#40)
### new features:
- check for new reference data on startup & update if necessary (#4)
- add additional line to IPD files to track the typeLoader version (#45)
### Improved user friendliness:
- Navigation area: doubleclick now opens SampleView or ProjectView (#51)
- enable editing of external sample IDs (#33)
- prohibit spaces in project names (#47)
### Adjustments for external users:
- enable and document running on Linux (#35)
- check whether ENA/IPD settings are filled out & create warning if not (#50)
- if blastn not found, generate helpful error message (#46)

## 2.0.1 (2018-10-12)
- Bugfix: if root\_path/_general could not be created during setup, create it during startup (#39)
- Bugfix: catch non-existing or empty config files during startup & exit gracefully (#36)
- improved user manual (#44)
- Bugfix: added table template .csv files to repo (#40)
- improved installer: (#43)
 - installer can now also be used to update an existing installation (this skips the creation of config files)
 - tidied-up headers etc.
 - now based on Python 3.6.6
 - cleaned uninstaller: delete previously leftover dirs (#41)

## 2.0.0 (2018-09-07) - official release to external!
