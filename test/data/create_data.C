void create_data() {

    TFile tf("data.root", "recreate");

    std::pair<double,double> range = {0, 100};
    int nbins = 100;

    // data set 1

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

    gRandom->SetSeed(0);
    h_data.FillRandom("f1", 200);
    h_data.FillRandom("f2", 150);
    h_data.FillRandom("f3", 300);

    f1.Write();
    f2.Write();
    f3.Write();
    h_f1.Write();
    h_f2.Write();
    h_f3.Write();
    h_f1_prior.Write();
    h_data.Write();

    // data set 2

    TF1 f4("f4", "gaus", range.first, range.second); f4.SetParameters(1, 50, 3);
    TF1 f5("f5", "[0]*exp([1]*x)", range.first, range.second); f5.SetParameters(1, -0.03);

    TH1D h_data_2("h_data_2", "Test Data 2", nbins, range.first, range.second);
    TH1D h_f4("h_f4", "Gauss 2", nbins, range.first, range.second);
    TH1D h_f5("h_f5", "Exp 1", nbins, range.first, range.second);

    for (int i = 1; i <= h_f4.GetNbinsX(); ++i) h_f4.SetBinContent(i, f4.Eval(h_f4.GetBinCenter(i)));
    for (int i = 1; i <= h_f5.GetNbinsX(); ++i) h_f5.SetBinContent(i, f5.Eval(h_f5.GetBinCenter(i)));

    h_data_2.FillRandom("f3", 300);
    h_data_2.FillRandom("f4", 200);
    h_data_2.FillRandom("f5", 5000);

    h_data_2.Write();
    h_f4.Write();
    h_f5.Write();

    // data set 3

    TF2 f6("f6", "gaus(x) + gaus(y)", range.first, range.second, range.first, range.second); f6.SetParameters(1, 50, 3);
    TF2 f7("f7", "[0]+[1]*(x+y)", range.first, range.second, range.first, range.second); f7.SetParameters(0, 10);

    TH2D h_data_3("h_data_3", "Test Data 3", nbins, range.first, range.second, nbins, range.first, range.second);
    TH2D h_f6("h_f6", "Gauss 3", nbins, range.first, range.second, nbins, range.first, range.second);
    TH2D h_f7("h_f7", "Bkg 2", nbins, range.first, range.second, nbins, range.first, range.second);

    for (int i = 1; i <= h_f6.GetNbinsX(); ++i) {
        for (int j = 1; j <= h_f6.GetNbinsY(); ++j) {
            h_f6.SetBinContent(i, j, f6.Eval(h_f6.GetXaxis()->GetBinCenter(i), h_f6.GetYaxis()->GetBinCenter(j)));
        }
    }
    for (int i = 1; i <= h_f7.GetNbinsX(); ++i) {
        for (int j = 1; j <= h_f7.GetNbinsY(); ++j) {
            h_f7.SetBinContent(i, j, f7.Eval(h_f7.GetXaxis()->GetBinCenter(i), h_f7.GetYaxis()->GetBinCenter(j)));
        }
    }

    h_data_3.FillRandom("f6", 500*500);
    h_data_3.FillRandom("f7", 1000*1000);

    h_data_3.Write();
    h_f6.Write();
    h_f7.Write();
}
