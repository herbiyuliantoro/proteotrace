# Releasing ProteoTrace

The repository is public. For a versioned release:

1. Set the same version in `pyproject.toml` and `proteotrace/__init__.py`.
2. Rebuild the figures and run the validation described in the repository
   documentation.
3. Review the changes, then create a GitHub release with a matching version
   tag.
4. If a DOI is needed, enable the repository in Zenodo and archive that
   release.

Add a paper citation to `CITATION.cff` after the paper or preprint has been
published and its DOI is available. Do not add a placeholder DOI.
