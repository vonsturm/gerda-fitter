void create_data() {

    std::pair<double,double> range = {0, 100};
    int nbins = 100;

    TF1 f1("f1", "gaus", range.first, range.second); f1.SetParameters(1, 20, 10);
    TF1 f2("f2", "gaus", range.first, range.second); f2.SetParameters(1, 80, 5);
    TF1 f3("f3", "[0]+[1]*x", range.first, range.second); f3.SetParameters(0, 10);

    TF1 f1_prior("f1_prior", "gaus", 0, 10); f1_prior.SetParameters(1, 0.005, 0.003);

    TH1D h_data("h_data", "Test Data", nbins, range.first, range.second);
    TH1D h_f1("h_f1", "Gauss 1", nbins, range.first, range.second);
    TH1D h_f2("h_f2", "Gauss 2", nbins, range.first, range.second);
    TH1D h_f3("h_f3", "Bkg 1", nbins, range.first, range.second);

    TH1D h_f1_prior("h_f1_prior", "Prior on Bkg 1", nbins, 0, 0.02);

    for (int i = 1; i <= h_f1.GetNbinsX(); ++i) h_f1.SetBinContent(i, f1.Eval(h_f1.GetBinCenter(i)));
    for (int i = 1; i <= h_f2.GetNbinsX(); ++i) h_f2.SetBinContent(i, f2.Eval(h_f2.GetBinCenter(i)));
    for (int i = 1; i <= h_f3.GetNbinsX(); ++i) h_f3.SetBinContent(i, f3.Eval(h_f3.GetBinCenter(i)));
    for (int i = 1; i <= h_f1_prior.GetNbinsX(); ++i) h_f1_prior.SetBinContent(i, f1_prior.Eval(h_f1_prior.GetBinCenter(i)));

    h_data.FillRandom("f1", 200);
    h_data.FillRandom("f2", 150);
    h_data.FillRandom("f3", 300);

    TFile tf("data.root", "recreate");
    f1.Write();
    f2.Write();
    f3.Write();
    h_f1.Write();
    h_f2.Write();
    h_f3.Write();
    h_f1_prior.Write();
    h_data.Write();
}
