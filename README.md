## simml-compare

This repository is a convenience to be able to compare versions 0.1.0, 0.2.0, and 0.3.0 of the core code for https://github.com/cran/simml 

## FAQs

### Why was this necessary when https://github.com/cran/simml contains tags for each version?

Because `simml-main.R` in version 0.1.0 was renamed `simml.main.R` in subsequent versions, but apparently _not_ using `git mv` to rename, so GitHub considers
them as completely different files, thus precluding comparison using the GitHub "diff" interface.
