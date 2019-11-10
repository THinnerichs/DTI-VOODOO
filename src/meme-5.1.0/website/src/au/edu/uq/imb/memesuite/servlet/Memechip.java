package au.edu.uq.imb.memesuite.servlet;

import au.edu.uq.imb.memesuite.data.*;
import au.edu.uq.imb.memesuite.db.MotifDB;
import au.edu.uq.imb.memesuite.db.MotifDBFile;
import au.edu.uq.imb.memesuite.servlet.util.*;
import au.edu.uq.imb.memesuite.template.HTMLSub;
import au.edu.uq.imb.memesuite.template.HTMLTemplate;
import au.edu.uq.imb.memesuite.util.FileCoord;
import au.edu.uq.imb.memesuite.util.JsonWr;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.activation.DataSource;
import javax.servlet.*;
import javax.servlet.http.*;

import static au.edu.uq.imb.memesuite.servlet.util.WebUtils.*;


public class Memechip extends SubmitJob<Memechip.Data> {
  private HTMLTemplate tmplMain;
  private HTMLTemplate tmplVerify;
  private ComponentHeader header;
  private ComponentMotifs motifs;
  private ComponentSequenceAlphabet alphabet;
  private ComponentSequences sequences;
  private ComponentSequences control;
  private ComponentBfile background;
  private ComponentJobDetails jobDetails;
  private ComponentAdvancedOptions universalOpts;
  private ComponentAdvancedOptions memeOpts;
  private ComponentAdvancedOptions dremeOpts;
  private ComponentAdvancedOptions centrimoOpts;
  private ComponentSubmitReset submitReset;
  private ComponentFooter footer;

  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web.memechip");
  
  protected class Data extends SubmitJob.JobData {
    public String email;
    public String description;
    public AlphabetDataSource alphabet;
    public SequenceDataSource sequences;
    public SequenceDataSource negSeq;
    public MotifInfo motifs;
    public String disc_mode;
    // universal options
    public boolean norc;
    public Background background;
    // MEME options
    String memeOptMode;
    int memeOptNMotifs;
    int memeOptMinW;
    int memeOptMaxW;
    Integer memeOptMinSites;
    Integer memeOptMaxSites;
    boolean memeOptPal; 	// look for palindromes
    boolean memeOptNorand;	// MEME won't randomize sequence order
    // DREME options
    double dremeOptE;
    Integer dremeOptM;
    // CentriMo options
    double centrimoOptScore;
    Integer centrimoOptMaxReg;
    double centrimoOptEThresh;
    boolean centrimoOptLocal;
    boolean centrimoOptStoreIds;

    @Override
    public void outputJson(JsonWr out) throws IOException {
      out.startObject();
      out.property("sequences", sequences);
      if (negSeq != null) out.property("negSeq", negSeq);
      out.property("disc_mode", disc_mode);
      out.property("motifs", motifs);
      out.property("norc", norc);
      out.property("background", background);
      out.property("memeOptMode", memeOptMode);
      out.property("memeOptNMotifs", memeOptNMotifs);
      out.property("memeOptMinW", memeOptMinW);
      out.property("memeOptMaxW", memeOptMaxW);
      if (memeOptMinSites != null) out.property("memeOptMinSites", memeOptMinSites);
      if (memeOptMaxSites != null) out.property("memeOptMaxSites", memeOptMaxSites);
      out.property("memeOptPal", memeOptPal);
      out.property("memeOptNorand", memeOptNorand);
      out.property("dremeOptE", dremeOptE);
      if (dremeOptM != null) out.property("dremeOptM", dremeOptM);
      out.property("centrimoOptScore", centrimoOptScore);
      if (centrimoOptMaxReg != null) out.property("centrimoOptMaxReg", centrimoOptMaxReg);
      out.property("centrimoOptEThresh", centrimoOptEThresh);
      out.property("centrimoOptLocal", centrimoOptLocal);
      out.property("centrimoOptStoreIds", centrimoOptStoreIds);
      out.endObject();
    }
    
    @Override
    public String email() {
      return email;
    }
  
    @Override
    public String description() {
      return description;  // generated code
    }

    @Override
    public String emailTemplate() {
      return tmplVerify.getSubtemplate("message").toString();
    }
  
    @Override
    public String cmd() {
      StringBuilder args = new StringBuilder();
      if (alphabet != null) {
        addArgs(args, "-alphf", alphabet.getName());
      } else {
        addArgs(args, "-alpha", sequences.guessAlphabet().name());
      }
      if (motifs instanceof MotifDataSource) {
        addArgs(args, "-upmotif", ((MotifDataSource) motifs).getName());
      }
      if (background.getSource() == Background.Source.FILE) {
        addArgs(args, "-bfile", background.getBfile().getName());
      } else {
        addArgs(args, "-order", background.getSource().getGeneratedOrder());
      }
      if (disc_mode == "psp") addArgs(args, "-psp-gen");
      if (negSeq != null) addArgs(args, "-neg", negSeq.getName());
      if (norc) addArgs(args, "-norc");
      // MEME specific arguments
      addArgs(args, "-meme-mod", memeOptMode);
      addArgs(args, "-meme-minw", memeOptMinW);
      addArgs(args, "-meme-maxw", memeOptMaxW);
      addArgs(args, "-meme-nmotifs", memeOptNMotifs);
      if (memeOptMinSites != null) addArgs(args, "-meme-minsites", memeOptMinSites);
      if (memeOptMaxSites != null) addArgs(args, "-meme-maxsites", memeOptMaxSites);
      if (memeOptPal) addArgs(args, "-meme-pal");
      if (memeOptNorand) addArgs(args, "-meme-norand");
      // DREME specific options
      addArgs(args, "-dreme-e", dremeOptE);
      if (dremeOptM != null) addArgs(args, "-dreme-m", dremeOptM);
      // CentriMo specific options
      if (centrimoOptLocal) addArgs(args, "-centrimo-local");
      addArgs(args, "-centrimo-score", centrimoOptScore);
      if (centrimoOptMaxReg != null) addArgs(args, "-centrimo-maxreg", centrimoOptMaxReg);
      addArgs(args, "-centrimo-ethresh", centrimoOptEThresh);
      if (!centrimoOptStoreIds) addArgs(args, "-centrimo-noseq");
      // sequences
      addArgs(args, sequences.getName());
      // motif databases
      if (motifs instanceof MotifDB) {
        for (MotifDBFile dbFile : ((MotifDB) motifs).getMotifFiles()) {
          addArgs(args, dbFile.getFileName());
        }
      }
      return args.toString();
    }
  
    @Override
    public List<DataSource> files() {
      ArrayList<DataSource> list = new ArrayList<DataSource>();
      if (alphabet != null) list.add(alphabet);
      list.add(sequences);
      if (negSeq != null) list.add(negSeq);
      if (motifs instanceof MotifDataSource) {
        list.add((MotifDataSource) motifs);
      }
      if (background.getSource() == Background.Source.FILE) {
        list.add(background.getBfile());
      }
      return list;
    }
  
    @Override
    public void cleanUp() {
      if (alphabet != null) {
        if (!alphabet.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              alphabet.getFile());
        }
      }
      if (sequences != null) {
        if (!sequences.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              sequences.getFile());
        }
      }
      if (negSeq != null) {
        if (!negSeq.getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              negSeq.getFile());
        }
      }
      if (motifs != null && motifs instanceof MotifDataSource) {
        if (!((MotifDataSource) motifs).getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              ((MotifDataSource) motifs).getFile());
        }
      }
      if (background.getSource() == Background.Source.FILE) {
        if (!background.getBfile().getFile().delete()) {
          logger.log(Level.WARNING, "Unable to delete temporary file " +
              background.getBfile().getFile());
        }
      }
    }
  }

  public Memechip() {
    super("MEMECHIP", "MEME-ChIP");
  }

  @Override
  public void init() throws ServletException {
    super.init();
    // load the templates
    tmplMain = cache.loadAndCache("/WEB-INF/templates/meme-chip.tmpl");
    tmplVerify = cache.loadAndCache("/WEB-INF/templates/meme-chip_verify.tmpl");
    header = new ComponentHeader(cache, msp.getVersion(), tmplMain.getSubtemplate("header"));
    alphabet = new ComponentSequenceAlphabet(context, tmplMain.getSubtemplate("alphabet"));
    motifs = new ComponentMotifs(context, tmplMain.getSubtemplate("motifs"));
    sequences = new ComponentSequences(context, tmplMain.getSubtemplate("sequences"));
    control = new ComponentSequences(context, tmplMain.getSubtemplate("control"));
    background = new ComponentBfile(context, tmplMain.getSubtemplate("bfile"));
    jobDetails = new ComponentJobDetails(cache);
    universalOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("universal_opts"));
    memeOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("meme_opts"));
    dremeOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("dreme_opts"));
    centrimoOpts = new ComponentAdvancedOptions(cache, tmplMain.getSubtemplate("centrimo_opts"));
    submitReset = new ComponentSubmitReset(cache, jobTable.getCount(), jobTable.getDuration());
    footer = new ComponentFooter(cache, msp);
  }

  @Override
  public String title() {
    return tmplVerify.getSubtemplate("title").toString();
  }

  @Override
  public String subtitle() {
    return tmplVerify.getSubtemplate("subtitle").toString();
  }

  @Override
  public String logoPath() {
    return tmplVerify.getSubtemplate("logo").toString();
  }

  @Override
  public String logoAltText() {
    return tmplVerify.getSubtemplate("alt").toString();
  }

  @Override
  protected void displayForm(HttpServletRequest request, HttpServletResponse response, long quotaMinWait) throws IOException {
    HTMLSub main = tmplMain.toSub();
    main.set("help", new HTMLSub[]{header.getHelp(), alphabet.getHelp(), motifs.getHelp(),
        sequences.getHelp(), background.getHelp(), jobDetails.getHelp(),
        universalOpts.getHelp(), submitReset.getHelp(), footer.getHelp()});
    main.set("header", header.getComponent());
    main.set("alphabet", alphabet.getComponent());
    main.set("motifs", motifs.getComponent());
    main.set("sequences", sequences.getComponent());
    main.set("control", control.getComponent());
    main.set("bfile", background.getComponent());
    main.set("job_details", jobDetails.getComponent());
    main.set("universal_opts", universalOpts.getComponent());
    main.set("meme_opts", memeOpts.getComponent());
    main.set("dreme_opts", dremeOpts.getComponent());
    main.set("centrimo_opts", centrimoOpts.getComponent());
    main.set("submit_reset", submitReset.getComponent(quotaMinWait));
    main.set("footer", footer.getComponent());
    response.setContentType("text/html; charset=UTF-8");
    main.output(response.getWriter());
  }

  @Override
  protected Data checkParameters(FeedbackHandler feedback,
      HttpServletRequest request) throws IOException, ServletException {
    // setup default file names
    FileCoord namer = new FileCoord();
    FileCoord.Name alphName = namer.createName("alphabet.alph");
    FileCoord.Name sequencesName = namer.createName("sequences.fa");
    FileCoord.Name negSeqName = namer.createName("control_sequences.fa");
    FileCoord.Name motifsName = namer.createName("motifs.meme");
    FileCoord.Name backgroundName = namer.createName("background");
    namer.createName("description");
    namer.createName("uuid");
    Alph alph = null;
    Data data = new Data();
    // get the job details
    data.email = jobDetails.getEmail(request, feedback);
    data.description = jobDetails.getDescription(request);
    // get the alphabet
    data.alphabet = alphabet.getAlphabet(alphName, request, feedback);
    if (data.alphabet != null) alph = data.alphabet.getAlph();
    // get the sequences
    data.sequences = (SequenceDataSource)sequences.getSequences(alph, sequencesName, request, feedback);
    if (alph == null && data.sequences != null) alph = data.sequences.guessAlphabet().getAlph();
    // get the discovery mode
    data.disc_mode = paramChoice(request, "disc_mode", "classic", "de", "psp");
    // get any negative sequences sequences
    if (data.disc_mode != "classic") {
      data.negSeq = (SequenceDataSource)control.getSequences(alph, negSeqName, request, feedback);
      // if (alph == null && data.negSeq != null) alph = data.negSeq.guessAlphabet().getAlph();
    }
    // get the motifs
    data.motifs =  motifs.getMotifs(alph, motifsName, request, feedback);
    // get the universal options
    data.norc = paramBool(request, "norc");
    data.background = background.getBfile(backgroundName, request, feedback);
    // get the MEME options
    data.memeOptMode = paramChoice(request, "meme_dist", "zoops", "oops", "anr");
    data.memeOptNMotifs = paramInteger(feedback, "MEME number of motifs", request, "meme_nmotifs", 0, null, 3);
    data.memeOptMinW = paramInteger(feedback, "MEME minimum motif width", request, "meme_minw", 2, 300, 6);
    data.memeOptMaxW = paramInteger(feedback, "MEME maximum motif width", request, "meme_maxw", 2, 300, 30);
    if (data.memeOptMinW > data.memeOptMaxW) {
      feedback.whine("The MEME minimum motif width is larger than the MEME maximum motif width.");
    }
    if (paramBool(request, "meme_minsites_enable")) {
      data.memeOptMinSites = paramInteger(feedback, "MEME minimum motif sites", request, "meme_minsites", 2, 600, 2);
    } else {
      data.memeOptMinSites = null;
    }
    if (paramBool(request, "meme_maxsites_enable")) {
      data.memeOptMaxSites = paramInteger(feedback, "MEME maximum motif sites", request, "meme_maxsites", 2, 600, 600);
      if (data.memeOptMinSites != null && data.memeOptMinSites > data.memeOptMaxSites) {
        feedback.whine("The MEME minimum motif sites is larger than the MEME maximum motif sites.");
      }
    } else {
      data.memeOptMaxSites = null;
    }
    data.memeOptPal = paramBool(request, "meme_pal");
    data.memeOptNorand = paramBool(request, "meme_norand");
    // get the DREME options
    data.dremeOptE = paramNumber(feedback, "DREME E-value threshold", request, "dreme_ethresh", 0.0, null, 0.05);
    if (paramBool(request, "dreme_nmotifs_enable")) {
      data.dremeOptM = paramInteger(feedback, "DREME motif count", request, "dreme_nmotifs", 0, null, 10);
    } else {
      data.dremeOptM = null;
    }
    // get the CentriMo options
    data.centrimoOptScore = paramNumber(feedback, "CentriMo match score threshold", request, "centrimo_score", null, null, 5.0);
    if (paramBool(request, "centrimo_maxreg_enable")) {
      data.centrimoOptMaxReg = paramInteger(feedback, "CentriMo maximum region", request, "centrimo_maxreg", 0, null, 200);
    } else {
      data.centrimoOptMaxReg = null;
    }
    data.centrimoOptEThresh = paramNumber(feedback, "CentriMo E-value threshold", request, "centrimo_ethresh", 0.0, null, 10.0);
    data.centrimoOptLocal = paramBool(request, "centrimo_local");
    data.centrimoOptStoreIds = paramBool(request, "centrimo_store_ids");
    return data;
  }

}

