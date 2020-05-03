void create_data() {

    std::pair<double,double> range = {0, 100};
    int nbins = 100;

    TF1 f1("f1", "gaus(0)", range.first, range.second); f1.SetParameters(1, 20, 10);
    TF1 f2("f2", "gaus(0)", range.first, range.second); f2.SetParameters(1, 80, 5);
    TF1 f3("f3", "[0]+[1]*x", range.first, range.second); f3.SetParameters(0, 10);

    TH1D h_data("h_data", "Test Data", nbins, range.first, range.second);
    TH1D h_f1("h_f1", "Gauss 1", nbins, range.first, range.second);
    TH1D h_f2("h_f2", "Gauss 2", nbins, range.first, range.second);
    TH1D h_f3("h_f3", "Bkg 1", nbins, range.first, range.second);

    for (int i = 1; i <= h_f1.GetNbinsX(); ++i) h_f1.SetBinContent(i, f1.Eval(h_f1.GetBinCenter(i)));
    for (int i = 1; i <= h_f2.GetNbinsX(); ++i) h_f2.SetBinContent(i, f2.Eval(h_f2.GetBinCenter(i)));
    for (int i = 1; i <= h_f3.GetNbinsX(); ++i) h_f3.SetBinContent(i, f3.Eval(h_f3.GetBinCenter(i)));

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
    h_data.Write();
}
