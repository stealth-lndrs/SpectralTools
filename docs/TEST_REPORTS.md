---
layout: default
title: "Test Reports"
---

# Test Reports

Use this page to track how to generate and retrieve the automated reports
produced by `scripts/run_full_report.jl`.

## Generating the Report

From the repository root:

```bash
julia --project=. --code-coverage=user scripts/run_full_report.jl
```

This prints a console summary **and** writes a markdown file under
`reports/test_report_<timestamp>.md`, e.g.
`reports/test_report_2025-12-12_161237.md`. The report includes:

- Pass/fail/error counts for every testset (including nested sets).
- Detailed listings of any failures/errors.
- Coverage per file, provided Julia was launched with `--code-coverage=user`.

## Viewing Reports Online

GitHub Pages can only serve files under `/docs`, so the `reports/` folder is not
directly exposed. Instead:

- Navigate to the repository on GitHub:
  `https://github.com/<user>/PoliSpectralTools.jl/tree/main/reports`
  (replace `<user>` with your account/organization).
- Open the most recent `test_report_<timestamp>.md` to view the rendered output
  in GitHub’s markdown viewer.

## Attaching to Releases/PRs

- After running the script, upload the latest report to your release notes or PR
  description so reviewers see the exact test/coverage snapshot.
- For extra visibility, paste the short console summary in the PR body and link
  to the markdown report (either via GitHub or by copying it into the PR).

## CI Integration

To automate report generation in GitHub Actions:

```yaml
steps:
  - uses: actions/checkout@v4
  - uses: julia-actions/setup-julia@v1
    with:
      version: '1.9'
  - run: julia --project=. -e 'using Pkg; Pkg.instantiate()'
  - run: julia --project=. --code-coverage=user scripts/run_full_report.jl
  - uses: actions/upload-artifact@v4
    with:
      name: test-report
      path: reports/test_report_*.md
```

This stores the markdown report as a CI artifact you can download from the
workflow summary page.
