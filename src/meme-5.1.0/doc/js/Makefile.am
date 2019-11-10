
SEDSPEC = -e 's%\@VERSION\@%$(VERSION)%g' \
          -e 's%\@CONTACT\@%$(CONTACT)%g' \
          -e 's%\@DEV_CONTACT\@%$(DEV_CONTACT)%g' \
          -e 's%\@SOURCE_URL\@%$(SOURCE_URL)%g' \
          -e 's%\@SITE_URL\@%$(URL)%g' \
          -e 's%\@ALT_URL\@%$(ALT_URL)%g' \
          -e 's%\@URL\@%$(URL)%g' \
          -e 's%\@PREV_URL\@%$(PREV_URL)%g' \
          -e 's%\@PREV_VER\@%$(PREV_VER)%g' \
          -e 's%\@NOTICES_URL\@%$(NOTICES_URL)%g' \
          -e 's%\@NEWS_URL\@%$(NEWS_URL)%g'


%.js : %.js.in
	$(SED) $(SEDSPEC) $< > $@

# files linked from the etc directory
SYMLINK_FILES = \
 constants.js.in \
 centrimo_graph.js \
 tgene_doc.js \
 delay_draw.js \
 eps_context.js \
 motif_logo.js \
 simple-shared-doc.js \
 spamo_graph.js \
 utilities.js

NORMAL_FILES = \
  Markdown.Converter.js \
  Markdown.Converter.LICENSE \
  manual.js \
  menu-path.js \
  shared-doc.js

# Included in distribution, but not installed.
IN_FILES = \
  constants.js.in \
  menu.js.in \
  menu-configure.js.in \
  menu-data.js.in

# Not included in distribution, but installed.
BUILT_FILES = \
  constants.js \
  menu.js \
  menu-configure.js \
  menu-data.js

BUILT_SOURCES = $(BUILT_FILES)

docjsdir = ${memedocdir}/js
docjs_DATA =  $(NORMAL_FILES) $(BUILT_FILES) $(SYMLINK_FILES)

CLEANFILES = $(BUILT_FILES) 
EXTRA_DIST = $(NORMAL_FILES) $(IN_FILES) $(SYMLINK_FILES)
