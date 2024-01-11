#####################################################
# 1. Properties for literature models
#####################################################

has_arxiv_id(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_id)
has_arxiv_doi(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_doi)
has_arxiv_link(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_link)
has_arxiv_model_equation_number(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_model_equation_number)
has_arxiv_model_page(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_model_page)
has_arxiv_model_section(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_model_section)
has_arxiv_version(m::AbstractFTheoryModel) = has_attribute(m, :arxiv_version)
has_associated_literature_models(m::AbstractFTheoryModel) = has_attribute(m, :associated_literature_models)
has_generating_sections(m::AbstractFTheoryModel) = has_attribute(m, :generating_sections)
has_journal_doi(m::AbstractFTheoryModel) = has_attribute(m, :journal_doi)
has_journal_link(m::AbstractFTheoryModel) = has_attribute(m, :journal_link)
has_journal_model_equation_number(m::AbstractFTheoryModel) = has_attribute(m, :journal_model_equation_number)
has_journal_model_page(m::AbstractFTheoryModel) = has_attribute(m, :journal_model_page)
has_journal_model_section(m::AbstractFTheoryModel) = has_attribute(m, :journal_model_section)
has_journal_pages(m::AbstractFTheoryModel) = has_attribute(m, :journal_pages)
has_journal_report_numbers(m::AbstractFTheoryModel) = has_attribute(m, :journal_report_numbers)
has_journal_volume(m::AbstractFTheoryModel) = has_attribute(m, :journal_volume)
has_journal_year(m::AbstractFTheoryModel) = has_attribute(m, :journal_year)
has_literature_identifier(m::AbstractFTheoryModel) = has_attribute(m, :literature_identifier)
has_model_description(m::AbstractFTheoryModel) = has_attribute(m, :model_description)
has_model_parameters(m::AbstractFTheoryModel) = has_attribute(m, :model_parameters)
has_paper_authors(m::AbstractFTheoryModel) = has_attribute(m, :paper_authors)
has_paper_buzzwords(m::AbstractFTheoryModel) = has_attribute(m, :paper_buzzwords)
has_paper_description(m::AbstractFTheoryModel) = has_attribute(m, :paper_description)
has_paper_title(m::AbstractFTheoryModel) = has_attribute(m, :paper_title)
has_related_literature_models(m::AbstractFTheoryModel) = has_attribute(m, :related_literature_models)
has_resolutions(m::AbstractFTheoryModel) = has_attribute(m, :resolutions)
has_resolution_generating_sections(m::AbstractFTheoryModel) = has_attribute(m, :resolution_generating_sections)
has_resolution_zero_sections(m::AbstractFTheoryModel) = has_attribute(m, :resolution_zero_sections)
has_weighted_resolutions(m::AbstractFTheoryModel) = has_attribute(m, :weighted_resolutions)
has_weighted_resolution_generating_sections(m::AbstractFTheoryModel) = has_attribute(m, :weighted_resolution_generating_sections)
has_weighted_resolution_zero_sections(m::AbstractFTheoryModel) = has_attribute(m, :weighted_resolution_zero_sections)
has_zero_section(m::AbstractFTheoryModel) = has_attribute(m, :zero_section)
