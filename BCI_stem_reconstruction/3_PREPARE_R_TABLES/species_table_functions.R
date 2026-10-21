# ========================================================================
# FUNCTIONS FOR THE SPECIES TABLE OF THE BCI R TABLES
# ========================================================================
# Sourced by 3_create_species_table_BCI.R and by tests/test_species_table.R,
# so the script and the tests run exactly the same code.
#   clean_text()          text without NA, outer or repeated spaces
#   full_name()           full name from its parts
#   latin_stem()          epithet without its gender ending
#   name_change_kind()    kind of change between two names
#   format_code_counts()  codes with their tree counts, "tremsp (84)"
#   absorbed_sentence()   note on a code no longer used, whose trees are
#                         now under another code
#   join_notes()          sentences of a note joined into one text
#   split_epithet()       epithet and infraspecific name of one text field
#   authority_key()       authority without its format (spaces, periods,
#                         capitals, accents), to compare authors only
#   name_parts_changed()  parts of a name that differ between two versions
# Every function is vectorised over its arguments (one element per code or
# per pair of codes). Requires data.table.
# ========================================================================

# Text without NA, outer spaces or repeated spaces, so that names are compared
# on their letters only ("" for NA)
clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("\\s+", " ", trimws(x))
}

# Full name from its parts (genus, epithet, infraspecific rank and epithet);
# empty parts are dropped, and a name without any part is ""
full_name <- function(...) clean_text(do.call(paste, lapply(list(...), clean_text)))

# Stem of a Latin epithet without its gender ending: -us/-a/-um, -is/-e,
# -er/-ra/-rum and -er/-era/-erum (adjectives), and -i/-ae and -ii/-iae
# (genitives of personal names, for a man or a woman). Two epithets that differ
# only in that ending share a stem: micrantha / micranthum, frondosum /
# frondosa, garwoodii / garwoodiae.
latin_stem <- function(x) sub("(erum|iae|era|rum|ii|ae|us|um|is|er|ra|a|e|i)$", "", clean_text(x))

# Kind of change between an old and a new name, given by their parts (infra =
# infraspecific rank and epithet, e.g. "var. grandiflora"). NA when the three
# parts are the same; otherwise the first of these that holds:
#   "morphospecies label"       the old name is not that of a species (Condit's
#                               IDLevel is not "species"), so its epithet is a
#                               label (sapisp: Sapium "broadleaf" -> "sp.")
#   "infraspecific name"        only the infraspecific name changed
#   "genus"                     the same epithet in another genus
#   "genus and epithet gender"  another genus, and the epithet changes only
#                               its gender ending (latin_stem): Capparidastrum
#                               frondosum -> Morisonia frondosa
#   "genus and epithet"         another genus and another epithet
#   "epithet gender"            the same genus, and the epithet changes only
#                               its gender ending: Trema micrantha ->
#                               micranthum (Trema is neuter)
#   "epithet spelling"          the same genus, and the epithet changes by at
#                               most max_letters letters otherwise (seibertii
#                               -> siebertii)
#   "epithet"                   the same genus, another epithet
# Letters = edit distance (adist). On 2026-10-15 the nearest different
# epithets of the same genus are 4 letters apart (crocata -> trigonata), so 2
# separates them. A change of the infraspecific name together with the genus
# or the epithet is labelled by the genus or the epithet. The labels describe
# the change only; they do not say which spelling is correct.
name_change_kind <- function(genus_old, epithet_old, infra_old,
                             genus_new, epithet_new, infra_new,
                             id_level_old = "species", max_letters = 2) {
  genus_old <- clean_text(genus_old)
  genus_new <- clean_text(genus_new)
  epithet_old <- clean_text(epithet_old)
  epithet_new <- clean_text(epithet_new)
  genus_changed <- genus_old != genus_new
  epithet_changed <- epithet_old != epithet_new
  infra_changed <- clean_text(infra_old) != clean_text(infra_new)
  epithet_gender <- epithet_changed & latin_stem(epithet_old) == latin_stem(epithet_new)
  letters_apart <- vapply(
    seq_along(epithet_old),
    function(k) drop(adist(epithet_old[k], epithet_new[k])),
    numeric(1)
  )
  id_level_old <- rep_len(as.character(id_level_old), length(genus_old))
  morphospecies <- !is.na(id_level_old) & id_level_old != "species"
  data.table::fcase(
    !genus_changed & !epithet_changed & !infra_changed, NA_character_,
    morphospecies, "morphospecies label",
    !genus_changed & !epithet_changed, "infraspecific name",
    genus_changed & !epithet_changed, "genus",
    genus_changed & epithet_gender, "genus and epithet gender",
    genus_changed, "genus and epithet",
    epithet_gender, "epithet gender",
    letters_apart <= max_letters, "epithet spelling",
    default = "epithet"
  )
}

# Codes with their tree counts, the largest group first (ties by code), e.g.
# "swars2 (26), protpa (2)"; NA when there are none. For one code at a time
# (call it by group).
format_code_counts <- function(codes, counts) {
  if (!length(codes)) {
    return(NA_character_)
  }
  o <- order(-counts, codes)
  paste0(codes[o], " (", counts[o], ")", collapse = ", ")
}

# Sentence for the note of a code that holds trees of a code no longer used:
# the old code and its name in the reference release, e.g. "Code tremsp,
# Trema unidentified in Condit et al. (2019), is no longer used." The number
# of trees is not repeated: it is in the column ReidentifiedFrom.
absorbed_sentence <- function(old_code, old_name, reference = "Condit et al. (2019)") {
  paste0("Code ", old_code, ", ", old_name, " in ", reference, ", is no longer used.")
}

# The sentences of a note joined into one text, row by row: each part is
# trimmed and ends with a full stop, empty parts are dropped, and a note
# without any part is NA. Every argument is a character vector with one
# element per row (NA or "" = no sentence).
join_notes <- function(...) {
  parts <- lapply(list(...), function(x) {
    x <- clean_text(x)
    data.table::fifelse(x == "" | grepl("[.!?]$", x), x, paste0(x, "."))
  })
  out <- clean_text(do.call(paste, parts))
  out[out == ""] <- NA_character_
  out
}

# Epithet and infraspecific name of a field that holds both, as the species
# column of the Flora de Panamá list does for the varieties ("simplex var.
# grandiflora", sometimes with two spaces): the first word is the epithet, the
# rest the infraspecific name with its rank ("" when there is none).
split_epithet <- function(x) {
  x <- clean_text(x)
  list(
    epithet = sub(" .*$", "", x),
    infra = data.table::fifelse(grepl(" ", x, fixed = TRUE), sub("^[^ ]+ ", "", x), "")
  )
}

# Authority reduced to its authors, so that two ways of writing the same
# authority compare equal: without spaces and periods, in lower case and
# without accents ("(L.) D. Don ex G. Don." and "(L.) D.Don ex G.Don";
# "(aubl.)" and "(Aubl.)"; "Müll. Arg." and "Mull.Arg."). Another author, an
# author added or dropped ("Engl." -> "(Kunth) Engl.") or another "ex" still
# differ. "" for NA.
authority_key <- function(x) {
  x <- tolower(gsub("[[:space:].]", "", clean_text(x)))
  chartr("áàâäãåéèêëíìîïóòôöõøúùûüñçý", "aaaaaaeeeeiiiioooooouuuuncy", x)
}

# Parts of a name that differ between version a and version b, as one text
# ("genus", "epithet", "infraspecific name", "authority", joined by ", "); NA
# when none differs. Genus, epithet and infraspecific name (with its rank)
# are compared letter by letter (clean_text: spaces aside); the authority by
# its authors only (authority_key).
name_parts_changed <- function(genus_a, epithet_a, infra_a, authority_a,
                               genus_b, epithet_b, infra_b, authority_b) {
  parts <- list(
    genus = clean_text(genus_a) != clean_text(genus_b),
    epithet = clean_text(epithet_a) != clean_text(epithet_b),
    `infraspecific name` = clean_text(infra_a) != clean_text(infra_b),
    authority = authority_key(authority_a) != authority_key(authority_b)
  )
  labels <- do.call(cbind, lapply(names(parts), function(p) data.table::fifelse(parts[[p]], p, NA_character_)))
  if (!length(labels)) {
    return(character(0))
  }
  out <- apply(labels, 1L, function(r) paste(r[!is.na(r)], collapse = ", "))
  out[out == ""] <- NA_character_
  out
}
