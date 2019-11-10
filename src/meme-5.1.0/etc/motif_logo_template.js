function motif_logo_template(inputs) {
  function _input(name) {
    if (typeof inputs[name] === "undefined") {
      throw new Error("Missing template variable: " + name);
    }
    return inputs[name];
  }
  return (
"%!PS-Adobe-3.0 EPSF-3.0\n" +
"%%Title: Sequence Logo : " + _input("TITLE") + "\n" +
"%%Creator: " + _input("CREATOR") + "\n" +
"%%CreationDate: " + _input("CREATIONDATE") + "\n" +
"%%BoundingBox:   0  0  " + _input("BOUNDINGWIDTH") + " " + _input("BOUNDINGHEIGHT") + " \n" +
"%%Pages: 0\n" +
"%%DocumentFonts: \n" +
"%%EndComments\n" +
"\n" +
"% ---- CONSTANTS ----\n" +
"\/cmfactor 72 2.54 div def % defines points -> cm conversion\n" +
"\/cm {cmfactor mul} bind def % defines centimeters\n" +
"\n" +
"% ---- VARIABLES ----\n" +
"\n" +
"% NA = Nucleic Acid, AA = Amino Acid\n" +
"\/logoType (" + _input("LOGOTYPE") + ") def \n" +
"\n" +
"\/logoTitle (" + _input("TITLE") + ") def\n" +
"\n" +
"% Dimensions in cm\n" +
"\/logoWidth " + _input("LOGOWIDTH") + " cm def\n" +
"\/logoHeight " + _input("LOGOLINEHEIGHT") + " cm def\n" +
"\/totalHeight " + _input("LOGOHEIGHT") + " cm def\n" +
"\n" +
"\/yaxis " + _input("YAXIS") + " def\n" +
"\/yaxisLabel (" + _input("YAXISLABEL") + ") def\n" +
"\/yaxisBits  " + _input("BARBITS") + " def % bits\n" +
"\/yaxisTicBits " + _input("TICBITS") + " def\n" +
"\n" +
"\/xaxis " + _input("NUMBERING") + " def\n" +
"\/xaxisLabel (" + _input("XAXISLABEL") + ") def\n" +
"\/showEnds (" + _input("SHOWENDS") + ") def \n" +
"\n" +
"\/showFineprint true def\n" +
"\/fineprint (" + _input("FINEPRINT") + ") def\n" +
"\n" +
"\/charsPerLine " + _input("CHARSPERLINE") + " def\n" +
"\n" +
"\/showingBox " + _input("SHOWINGBOX") + " def    \n" +
"\/shrinking false def   % true falses\n" +
"\/shrink  1.0 def\n" +
"\/outline " + _input("OUTLINE") + " def\n" +
"\n" +
"\/IbeamFraction  " + _input("ERRORBARFRACTION") + " def\n" +
"\/IbeamGray      0.50 def\n" +
"\/IbeamLineWidth 0.5 def\n" +
"\n" +
"\/fontsize       " + _input("FONTSIZE") + " def\n" +
"\/titleFontsize  " + _input("TITLEFONTSIZE") + " def\n" +
"\/smallFontsize  " + _input("SMALLFONTSIZE") + " def\n" +
"\n" +
"\/topMargin      " + _input("TOPMARGIN") + " cm def\n" +
"\/bottomMargin   " + _input("BOTTOMMARGIN") + " cm def\n" +
"\n" +
"\/defaultColor [0 0 0] def \n" +
"\n" +
_input("COLORDICT") + "\n" +
"\n" +
"\/colorDict fullColourDict def\n" +
"\n" +
"% ---- DERIVED PARAMETERS ----\n" +
"\n" +
"\/leftMargin\n" +
"  fontsize 3.5 mul\n" +
"\n" +
"def \n" +
"\n" +
"\/rightMargin \n" +
"  %Add extra room if showing ends\n" +
"  showEnds (false) eq { fontsize}{fontsize 1.5 mul} ifelse\n" +
"def\n" +
"\n" +
"\/yaxisHeight \n" +
"  logoHeight \n" +
"  bottomMargin sub  \n" +
"  topMargin sub\n" +
"def\n" +
"\n" +
"\/ticWidth fontsize 2 div def\n" +
"\n" +
"\/pointsPerBit yaxisHeight yaxisBits div  def\n" +
"\n" +
"\/stackMargin 1 def\n" +
"\n" +
"% Do not add space aroung characters if characters are boxed\n" +
"\/charRightMargin \n" +
"  showingBox { 0.0 } {stackMargin} ifelse\n" +
"def\n" +
"\n" +
"\/charTopMargin \n" +
"  showingBox { 0.0 } {stackMargin} ifelse\n" +
"def\n" +
"\n" +
"\/charWidth\n" +
"  logoWidth\n" +
"  leftMargin sub\n" +
"  rightMargin sub\n" +
"  charsPerLine div\n" +
"  charRightMargin sub\n" +
"def\n" +
"\n" +
"\/charWidth4 charWidth 4 div def\n" +
"\/charWidth2 charWidth 2 div def\n" +
"\n" +
"\/stackWidth \n" +
"  charWidth charRightMargin add\n" +
"def\n" +
" \n" +
"\/numberFontsize \n" +
"  fontsize charWidth lt {fontsize}{charWidth} ifelse\n" +
"def\n" +
"\n" +
"% movements to place 5'\/N and 3'\/C symbols\n" +
"\/leftEndDeltaX  fontsize neg         def\n" +
"\/leftEndDeltaY  fontsize 1.5 mul neg def\n" +
"\/rightEndDeltaX fontsize 0.25 mul     def\n" +
"\/rightEndDeltaY leftEndDeltaY        def\n" +
"\n" +
"% Outline width is proporional to charWidth, \n" +
"% but no less that 1 point\n" +
"\/outlinewidth \n" +
"  charWidth 32 div dup 1 gt  {}{pop 1} ifelse\n" +
"def\n" +
"\n" +
"\n" +
"% ---- PROCEDURES ----\n" +
"\n" +
"\/StartLogo { \n" +
"  % Save state\n" +
"  save \n" +
"  gsave \n" +
"\n" +
"  % Print Logo Title, top center \n" +
"  gsave \n" +
"    SetStringFont\n" +
"\n" +
"    logoWidth 2 div\n" +
"    logoTitle\n" +
"    stringwidth pop 2 div sub\n" +
"    totalHeight\n" +
"    titleFontsize sub\n" +
"    moveto\n" +
"\n" +
"    logoTitle\n" +
"    show\n" +
"  grestore\n" +
"\n" +
"  % Print X-axis label, bottom center\n" +
"  gsave\n" +
"    SetStringFont\n" +
"\n" +
"    logoWidth 2 div\n" +
"    xaxisLabel\n" +
"    stringwidth pop 2 div sub\n" +
"    0\n" +
"    titleFontsize 3 div\n" +
"    add\n" +
"    moveto\n" +
"\n" +
"    xaxisLabel\n" +
"    show\n" +
"  grestore\n" +
"\n" +
"  % Show Fine Print\n" +
"  showFineprint {\n" +
"    gsave\n" +
"      SetSmallFont\n" +
"      logoWidth\n" +
"        fineprint stringwidth pop sub\n" +
"        smallFontsize sub\n" +
"          smallFontsize 3 div\n" +
"      moveto\n" +
"    \n" +
"      fineprint show\n" +
"    grestore\n" +
"  } if\n" +
"\n" +
"  % Move to lower left corner of last line, first stack\n" +
"  leftMargin bottomMargin translate\n" +
"\n" +
"  % Move above first line ready for StartLine \n" +
"  0 totalHeight translate\n" +
"\n" +
"  SetLogoFont\n" +
"} bind def\n" +
"\n" +
"\/EndLogo { \n" +
"  grestore \n" +
"  showpage \n" +
"  restore \n" +
"} bind def\n" +
"\n" +
"\n" +
"\/StartLine { \n" +
"  % move down to the bottom of the line:\n" +
"  0 logoHeight neg translate\n" +
"  \n" +
"  gsave \n" +
"    yaxis { MakeYaxis } if\n" +
"    xaxis { showEnds (true) eq {ShowLeftEnd} if } if\n" +
"} bind def\n" +
"\n" +
"\/EndLine{ \n" +
"    xaxis { showEnds (true) eq {ShowRightEnd} if } if\n" +
"  grestore \n" +
"} bind def\n" +
"\n" +
"\n" +
"\/MakeYaxis {\n" +
"  gsave    \n" +
"    stackMargin neg 0 translate\n" +
"    ShowYaxisBar\n" +
"    ShowYaxisLabel\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowYaxisBar { \n" +
"  gsave  \n" +
"    SetStringFont\n" +
"\n" +
"    \/str 10 string def % string to hold number  \n" +
"    \/smallgap stackMargin 2 div def\n" +
"\n" +
"    % Draw first tic and bar\n" +
"    gsave    \n" +
"      ticWidth neg 0 moveto \n" +
"      ticWidth 0 rlineto \n" +
"      0 yaxisHeight rlineto\n" +
"      stroke\n" +
"    grestore\n" +
"\n" +
"   \n" +
"    % Draw the tics\n" +
"    % initial increment limit proc for\n" +
"    0 yaxisTicBits yaxisBits abs %cvi\n" +
"    {\/loopnumber exch def\n" +
"\n" +
"      % convert the number coming from the loop to a string\n" +
"      % and find its width\n" +
"      loopnumber 10 str cvrs\n" +
"      \/stringnumber exch def % string representing the number\n" +
"\n" +
"      stringnumber stringwidth pop\n" +
"      \/numberwidth exch def % width of number to show\n" +
"\n" +
"      \/halfnumberheight\n" +
"         stringnumber CharBoxHeight 2 div\n" +
"      def\n" +
"\n" +
"      numberwidth % move back width of number\n" +
"      neg loopnumber pointsPerBit mul % shift on y axis\n" +
"      halfnumberheight sub % down half the digit\n" +
"\n" +
"      moveto % move back the width of the string\n" +
"\n" +
"      ticWidth neg smallgap sub % Move back a bit more  \n" +
"      0 rmoveto % move back the width of the tic  \n" +
"\n" +
"      stringnumber show\n" +
"      smallgap 0 rmoveto % Make a small gap  \n" +
"\n" +
"      % now show the tic mark\n" +
"      0 halfnumberheight rmoveto % shift up again\n" +
"      ticWidth 0 rlineto\n" +
"      stroke\n" +
"    } for\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\/ShowYaxisLabel {\n" +
"  gsave\n" +
"    SetStringFont\n" +
"\n" +
"    % How far we move left depends on the size of\n" +
"    % the tic labels.\n" +
"    \/str 10 string def % string to hold number  \n" +
"    yaxisBits yaxisTicBits div cvi yaxisTicBits mul \n" +
"    str cvs stringwidth pop\n" +
"    ticWidth 1.5 mul  add neg  \n" +
"\n" +
"\n" +
"    yaxisHeight\n" +
"    yaxisLabel stringwidth pop\n" +
"    sub 2 div\n" +
"\n" +
"    translate\n" +
"    90 rotate\n" +
"    0 0 moveto\n" +
"    yaxisLabel show\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/StartStack {  % <stackNumber> startstack\n" +
"  xaxis {MakeNumber}{pop} ifelse\n" +
"  gsave\n" +
"} bind def\n" +
"\n" +
"\/EndStack {\n" +
"  grestore\n" +
"  stackWidth 0 translate\n" +
"} bind def\n" +
"\n" +
"\n" +
"% Draw a character whose height is proportional to symbol bits\n" +
"\/MakeSymbol{ % charbits character MakeSymbol\n" +
"  gsave\n" +
"    \/char exch def\n" +
"    \/bits exch def\n" +
"\n" +
"    \/bitsHeight \n" +
"       bits pointsPerBit mul \n" +
"    def\n" +
"\n" +
"    \/charHeight \n" +
"       bitsHeight charTopMargin sub\n" +
"       dup \n" +
"       0.0 gt {}{pop 0.0} ifelse % if neg replace with zero \n" +
"    def \n" +
" \n" +
"    charHeight 0.0 gt {\n" +
"      char SetColor\n" +
"      charWidth charHeight char ShowChar\n" +
"\n" +
"      showingBox { % Unfilled box\n" +
"        0 0 charWidth charHeight false ShowBox\n" +
"      } if\n" +
"\n" +
"\n" +
"    } if\n" +
"\n" +
"  grestore\n" +
"\n" +
"  0 bitsHeight translate \n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowChar { % <width> <height> <char> ShowChar\n" +
"  gsave\n" +
"    \/tc exch def    % The character\n" +
"    \/ysize exch def % the y size of the character\n" +
"    \/xsize exch def % the x size of the character\n" +
"\n" +
"    \/xmulfactor 1 def \n" +
"    \/ymulfactor 1 def\n" +
"    \/limmulfactor 0.01 def\n" +
"    \/drawable true def\n" +
"\n" +
"  \n" +
"    % if ysize is negative, make everything upside down!\n" +
"    ysize 0 lt {\n" +
"      % put ysize normal in this orientation\n" +
"      \/ysize ysize abs def\n" +
"      xsize ysize translate\n" +
"      180 rotate\n" +
"    } if\n" +
"\n" +
"    shrinking {\n" +
"      xsize 1 shrink sub 2 div mul\n" +
"        ysize 1 shrink sub 2 div mul translate \n" +
"\n" +
"      shrink shrink scale\n" +
"    } if\n" +
"\n" +
"    % Calculate the font scaling factors\n" +
"    % Loop twice to catch small correction due to first scaling\n" +
"    2 {\n" +
"      gsave\n" +
"        xmulfactor ymulfactor scale\n" +
"      \n" +
"        ysize % desired size of character in points\n" +
"        tc CharBoxHeight \n" +
"        dup 0.0 ne {\n" +
"          div % factor by which to scale up the character\n" +
"          \/ymulfactor exch def\n" +
"        } % end if\n" +
"        {pop pop}\n" +
"        ifelse\n" +
"\n" +
"        xsize % desired size of character in points\n" +
"        tc CharBoxWidth  \n" +
"        dup 0.0 ne {\n" +
"          div % factor by which to scale up the character\n" +
"          \/xmulfactor exch def\n" +
"        } % end if\n" +
"        {pop pop}\n" +
"        ifelse\n" +
"      grestore\n" +
"      % if the multiplication factors get too small we need to avoid a crash\n" +
"      xmulfactor limmulfactor lt {\n" +
"        \/xmulfactor 1 def\n" +
"        \/drawable false def\n" +
"      } if\n" +
"      ymulfactor limmulfactor lt {\n" +
"        \/ymulfactor 1 def\n" +
"        \/drawable false def\n" +
"      } if\n" +
"    } repeat\n" +
"\n" +
"    % Adjust horizontal position if the symbol is an I\n" +
"    tc (I) eq {\n" +
"      charWidth 2 div % half of requested character width\n" +
"      tc CharBoxWidth 2 div % half of the actual character\n" +
"      sub 0 translate\n" +
"      % Avoid x scaling for I \n" +
"      \/xmulfactor 1 def \n" +
"    } if\n" +
"\n" +
"\n" +
"    % ---- Finally, draw the character\n" +
"    drawable { \n" +
"      newpath\n" +
"      xmulfactor ymulfactor scale\n" +
"\n" +
"      % Move lower left corner of character to start point\n" +
"      tc CharBox pop pop % llx lly : Lower left corner\n" +
"      exch neg exch neg\n" +
"      moveto\n" +
"\n" +
"      outline {  % outline characters:\n" +
"        outlinewidth setlinewidth\n" +
"        tc true charpath\n" +
"        gsave 1 setgray fill grestore\n" +
"        clip stroke\n" +
"      } { % regular characters\n" +
"        tc show\n" +
"      } ifelse\n" +
"    } if\n" +
"\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowBox { % x1 y1 x2 y2 filled ShowBox\n" +
"  gsave\n" +
"    \/filled exch def \n" +
"    \/y2 exch def\n" +
"    \/x2 exch def\n" +
"    \/y1 exch def\n" +
"    \/x1 exch def\n" +
"    newpath\n" +
"    x1 y1 moveto\n" +
"    x2 y1 lineto\n" +
"    x2 y2 lineto\n" +
"    x1 y2 lineto\n" +
"    closepath\n" +
"\n" +
"    clip\n" +
"    \n" +
"    filled {\n" +
"      fill\n" +
"    }{ \n" +
"      0 setgray stroke   \n" +
"    } ifelse\n" +
"\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/MakeNumber { % number MakeNumber\n" +
"  gsave\n" +
"    SetNumberFont\n" +
"    stackWidth 0 translate\n" +
"    90 rotate % rotate so the number fits\n" +
"    dup stringwidth pop % find the length of the number\n" +
"    neg % prepare for move\n" +
"    stackMargin sub % Move back a bit\n" +
"    charWidth (0) CharBoxHeight % height of numbers\n" +
"    sub 2 div %\n" +
"    moveto % move back to provide space\n" +
"    show\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/Ibeam{ % heightInBits Ibeam\n" +
"  gsave\n" +
"    % Make an Ibeam of twice the given height in bits\n" +
"    \/height exch  pointsPerBit mul def \n" +
"    \/heightDRAW height IbeamFraction mul def\n" +
"\n" +
"    IbeamLineWidth setlinewidth\n" +
"    IbeamGray setgray \n" +
"\n" +
"    charWidth2 height neg translate\n" +
"    ShowIbar\n" +
"    newpath\n" +
"      0 0 moveto\n" +
"      0 heightDRAW rlineto\n" +
"    stroke\n" +
"    newpath\n" +
"      0 height moveto\n" +
"      0 height rmoveto\n" +
"      currentpoint translate\n" +
"    ShowIbar\n" +
"    newpath\n" +
"    0 0 moveto\n" +
"    0 heightDRAW neg rlineto\n" +
"    currentpoint translate\n" +
"    stroke\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowIbar { % make a horizontal bar\n" +
"  gsave\n" +
"    newpath\n" +
"      charWidth4 neg 0 moveto\n" +
"      charWidth4 0 lineto\n" +
"    stroke\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowLeftEnd {\n" +
"  gsave\n" +
"    SetStringFont\n" +
"    leftEndDeltaX leftEndDeltaY moveto\n" +
"    logoType (NA) eq {(5) show ShowPrime} if\n" +
"    logoType (AA) eq {(N) show} if\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowRightEnd { \n" +
"  gsave\n" +
"    SetStringFont\n" +
"    rightEndDeltaX rightEndDeltaY moveto\n" +
"    logoType (NA) eq {(3) show ShowPrime} if\n" +
"    logoType (AA) eq {(C) show} if\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"\/ShowPrime {\n" +
"  gsave\n" +
"    SetPrimeFont\n" +
"    (\\242) show \n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
" \n" +
"\/SetColor{ % <char> SetColor\n" +
"  dup colorDict exch known {\n" +
"    colorDict exch get aload pop setrgbcolor\n" +
"  } {\n" +
"    pop\n" +
"    defaultColor aload pop setrgbcolor\n" +
"  } ifelse \n" +
"} bind def\n" +
"\n" +
"% define fonts\n" +
"\/SetTitleFont {\/Times-Bold findfont titleFontsize scalefont setfont} bind def\n" +
"\/SetLogoFont  {\/Helvetica-Bold findfont charWidth  scalefont setfont} bind def\n" +
"\/SetStringFont{\/Helvetica-Bold findfont fontsize scalefont setfont} bind def\n" +
"\/SetPrimeFont {\/Symbol findfont fontsize scalefont setfont} bind def\n" +
"\/SetSmallFont {\/Helvetica findfont smallFontsize scalefont setfont} bind def\n" +
"\n" +
"\/SetNumberFont {\n" +
"    \/Helvetica-Bold findfont \n" +
"    numberFontsize\n" +
"    scalefont\n" +
"    setfont\n" +
"} bind def\n" +
"\n" +
"%Take a single character and return the bounding box\n" +
"\/CharBox { % <char> CharBox <lx> <ly> <ux> <uy>\n" +
"  gsave\n" +
"    newpath\n" +
"    0 0 moveto\n" +
"    % take the character off the stack and use it here:\n" +
"    true charpath \n" +
"    flattenpath \n" +
"    pathbbox % compute bounding box of 1 pt. char => lx ly ux uy\n" +
"    % the path is here, but toss it away ...\n" +
"  grestore\n" +
"} bind def\n" +
"\n" +
"\n" +
"% The height of a characters bounding box\n" +
"\/CharBoxHeight { % <char> CharBoxHeight <num>\n" +
"  CharBox\n" +
"  exch pop sub neg exch pop\n" +
"} bind def\n" +
"\n" +
"\n" +
"% The width of a characters bounding box\n" +
"\/CharBoxWidth { % <char> CharBoxHeight <num>\n" +
"  CharBox\n" +
"  pop exch pop sub neg \n" +
"} bind def\n" +
"\n" +
"% Set the colour scheme to be faded to indicate trimming\n" +
"\/MuteColour {\n" +
"  \/colorDict mutedColourDict def\n" +
"} def\n" +
"\n" +
"% Restore the colour scheme to the normal colours\n" +
"\/RestoreColour {\n" +
"  \/colorDict fullColourDict def\n" +
"} def\n" +
"\n" +
"% Draw the background for a trimmed section\n" +
"% takes the number of columns as a parameter\n" +
"\/DrawTrimBg { % <num> DrawTrimBox\n" +
"  \/col exch def\n" +
"  \n" +
"  \/boxwidth \n" +
"    col stackWidth mul \n" +
"  def\n" +
" \n" +
"  gsave\n" +
"    0.97 setgray\n" +
"\n" +
"    newpath\n" +
"    0 0 moveto\n" +
"    boxwidth 0 rlineto\n" +
"    0 yaxisHeight rlineto\n" +
"    0 yaxisHeight lineto\n" +
"    closepath\n" +
"    \n" +
"    fill\n" +
"  grestore\n" +
"} def\n" +
"\n" +
"\/DrawTrimEdge {\n" +
"  gsave\n" +
"    0.2 setgray\n" +
"    [2] 0 setdash\n" +
"\n" +
"    newpath\n" +
"    0 0 moveto\n" +
"    0 yaxisHeight lineto\n" +
"    \n" +
"    stroke\n" +
"\n" +
"} def\n" +
"\n" +
"\n" +
"% Deprecated names\n" +
"\/startstack {StartStack} bind  def\n" +
"\/endstack {EndStack}     bind def\n" +
"\/makenumber {MakeNumber} bind def\n" +
"\/numchar { MakeSymbol }  bind def\n" +
"\n" +
"%%EndProlog\n" +
"\n" +
"%%Page: 1 1\n" +
"StartLogo\n" +
"\n" +
_input("DATA") + "\n" +
"\n" +
"EndLogo\n" +
"\n" +
"%%EOF\n"
  );
}