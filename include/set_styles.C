void gemLines(TLine *gem_line) {

	gem_line->SetLineStyle(10);
    gem_line->SetLineColor(6);
    gem_line->SetLineWidth(2);
    gem_line->Draw("same");

    return 0;
}