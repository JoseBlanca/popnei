//! The statistics of the variants and of the individuals, per population.
//!
//! A population is a named set of individuals that a calculation treats as
//! a group, and every statistic of this module is calculated for each
//! population over its individuals alone. [`Pops`] is what a pass works
//! with: the name of each population and the indices of its individuals
//! among the individuals the pass gives, which are those of the source
//! after the filter of individuals of `docs/specs/filters.md` when the
//! variants carry one.
//!
//! `docs/specs/stats.md` has the design, and the row `stats` of section 9
//! of `docs/architecture.md` where the module sits.

use crate::error::{Error, Result};
use crate::filters::resolve_individuals;

/// The name of the one population of a calculation that was given no
/// populations, inherited from pyNei's `DEF_POP_NAME`.
pub const DEFAULT_POP_NAME: &str = "pop";

/// One population: its name and the indices of its individuals.
#[derive(Debug)]
struct Pop {
    name: String,
    /// The index of each individual among those of the reader, in the order
    /// the user named them.
    individuals: Vec<usize>,
    /// Whether those indices are every individual of the reader, in its
    /// order.
    is_all: bool,
}

/// The populations a pass calculates its statistics for: each one's name
/// and the indices of its individuals among the individuals of the reader.
///
/// The populations are in the order they were given, which is the order the
/// keys of a user's `pops` iterate in, and every result of this module
/// holds its values in that order. An individual can be in two
/// populations, and is in each of them once.
#[derive(Debug)]
pub struct Pops {
    pops: Vec<Pop>,
}

impl Pops {
    /// One population, named [`DEFAULT_POP_NAME`], of every individual of a
    /// reader of `num_individuals`, in the order of the reader: what a
    /// calculation works over when the user named no population.
    #[must_use]
    pub fn all(num_individuals: usize) -> Pops {
        Pops {
            pops: vec![Pop {
                name: DEFAULT_POP_NAME.to_owned(),
                individuals: (0..num_individuals).collect(),
                is_all: true,
            }],
        }
    }

    /// The populations a user named, each name looked up among
    /// `individuals`, the individuals the pass gives.
    ///
    /// `pops` is the name of each population with the names of its
    /// individuals, in the order the user gave them; it comes from a dict
    /// in Python and from an object in TypeScript, which hold each
    /// population name once.
    ///
    /// # Errors
    ///
    /// A name that is not one of `individuals`, a name that is twice in one
    /// population, a population that names no individual, and no population
    /// at all. The first three name the population, and the first two the
    /// name the user wrote.
    pub fn from_names(pops: &[(String, Vec<String>)], individuals: &[String]) -> Result<Pops> {
        if pops.is_empty() {
            return Err(Error::NoPop);
        }
        let mut of_the_pass = Vec::with_capacity(pops.len());
        for (name, named) in pops {
            let of_the_pop =
                resolve_individuals(named, individuals).map_err(|error| in_the_pop(error, name))?;
            of_the_pass.push(Pop {
                name: name.clone(),
                is_all: every_individual_in_order(&of_the_pop, individuals.len()),
                individuals: of_the_pop,
            });
        }
        Ok(Pops { pops: of_the_pass })
    }

    /// How many populations there are, which is 1 when the user named none.
    #[must_use]
    pub fn len(&self) -> usize {
        self.pops.len()
    }

    /// Whether there is no population, which the constructors do not build:
    /// `pops` with no population is refused, and [`Pops::all`] gives one.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pops.is_empty()
    }

    /// The name of one population, as the user wrote it.
    ///
    /// `pop` is a population of `0..len()`, which is how a caller walks
    /// them; a number at or beyond `len()` is no population of this and has
    /// no name here.
    #[must_use]
    pub fn name(&self, pop: usize) -> &str {
        self.pops.get(pop).map_or("", |pop| pop.name.as_str())
    }

    /// The index of each individual of one population among the individuals
    /// of the reader, in the order the user named them.
    ///
    /// `pop` is a population of `0..len()`; a number at or beyond `len()`
    /// is no population of this and has no individual here.
    #[must_use]
    pub fn individuals(&self, pop: usize) -> &[usize] {
        self.pops
            .get(pop)
            .map_or(&[][..], |pop| pop.individuals.as_slice())
    }

    /// Whether the population is every individual of the reader in the
    /// order of the reader, which a caller that reads a row of genotypes as
    /// it is asks before it counts.
    ///
    /// `pop` is a population of `0..len()`; a number at or beyond `len()`
    /// is no population of this and is not every individual.
    #[must_use]
    pub fn is_all(&self, pop: usize) -> bool {
        self.pops.get(pop).is_some_and(|pop| pop.is_all)
    }
}

/// The error of [`resolve_individuals`] over the names of one population,
/// with the population they were given in.
///
/// The filter of individuals and a population are given the same three
/// refusals by the same function, and what a user has to look at differs:
/// there the names are the whole of what they wrote, and here they are the
/// names of one of their populations, which is the one the message has to
/// name. An error of another kind, which that function does not give,
/// travels on as it is.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "`resolve_individuals` fails with these three cases alone, and an error of \
              any other case of the crate is left as it is rather than named after a \
              population it may have nothing to do with"
)]
fn in_the_pop(error: Error, pop: &str) -> Error {
    match error {
        Error::IndividualNotInTheSource { name } => Error::IndividualOfAPopNotInThePass {
            pop: pop.to_owned(),
            name,
        },
        Error::IndividualNamedTwice { name } => Error::IndividualNamedTwiceInAPop {
            pop: pop.to_owned(),
            name,
        },
        Error::NoIndividualNamed => Error::PopWithNoIndividual {
            pop: pop.to_owned(),
        },
        of_another_kind => of_another_kind,
    }
}

/// Whether `individuals` is every individual of a reader of
/// `num_individuals`, in the order of the reader: the population a caller
/// counts by reading a row of genotypes as it is.
fn every_individual_in_order(individuals: &[usize], num_individuals: usize) -> bool {
    individuals.len() == num_individuals
        && individuals
            .iter()
            .enumerate()
            .all(|(of_the_reader, individual)| of_the_reader == *individual)
}

#[cfg(test)]
mod pops {
    use super::{DEFAULT_POP_NAME, Pops};
    use crate::error::Error;

    /// The five diploid individuals of the worked example of
    /// `docs/specs/filters.md`, in the order of the source.
    fn the_five_individuals() -> Vec<String> {
        ["i1", "i2", "i3", "i4", "i5"]
            .iter()
            .map(|name| (*name).to_owned())
            .collect()
    }

    /// One population as `from_names` takes it, written as a user writes
    /// it.
    fn pop_of(name: &str, individuals: &[&str]) -> (String, Vec<String>) {
        (
            name.to_owned(),
            individuals
                .iter()
                .map(|individual| (*individual).to_owned())
                .collect(),
        )
    }

    /// The two populations of the worked example of "How it is verified"
    /// of the per variant distributions, pop1 of i1 and i2 and pop2 of i3,
    /// i4 and i5, which are the individuals 0 and 1 and the individuals 2,
    /// 3 and 4 of the source. Neither is every individual of the source.
    #[test]
    fn the_indices_of_the_two_populations_of_the_worked_example() {
        let pops = Pops::from_names(
            &[
                pop_of("pop1", &["i1", "i2"]),
                pop_of("pop2", &["i3", "i4", "i5"]),
            ],
            &the_five_individuals(),
        )
        .unwrap();
        assert_eq!(pops.len(), 2);
        assert!(!pops.is_empty());
        assert_eq!(pops.name(0), "pop1");
        assert_eq!(pops.individuals(0), [0, 1]);
        assert!(!pops.is_all(0));
        assert_eq!(pops.name(1), "pop2");
        assert_eq!(pops.individuals(1), [2, 3, 4]);
        assert!(!pops.is_all(1));
    }

    /// The order of the populations is the one they were given in, and the
    /// order of the individuals of each the one the user named them in,
    /// which is not the order of the source here.
    #[test]
    fn the_populations_and_their_individuals_are_in_the_order_they_were_given() {
        let pops = Pops::from_names(
            &[pop_of("second", &["i5", "i1"]), pop_of("first", &["i3"])],
            &the_five_individuals(),
        )
        .unwrap();
        assert_eq!(pops.name(0), "second");
        assert_eq!(pops.individuals(0), [4, 0]);
        assert_eq!(pops.name(1), "first");
        assert_eq!(pops.individuals(1), [2]);
    }

    /// The name a user wrote and the population it is in are what they
    /// have to look at, so the error carries both.
    #[test]
    fn a_name_that_is_not_an_individual_is_refused_with_its_population() {
        let error = Pops::from_names(&[pop_of("pop1", &["i1", "nope"])], &the_five_individuals())
            .unwrap_err();
        assert!(
            matches!(&error, Error::IndividualOfAPopNotInThePass { pop, name }
                if pop == "pop1" && name == "nope"),
            "{error:?}"
        );
    }

    /// pyNei counts an individual named twice twice, and popnei refuses it:
    /// a count that is wrong and says nothing is what popnei never gives.
    #[test]
    fn a_name_twice_in_one_population_is_refused_with_its_population() {
        let error = Pops::from_names(
            &[pop_of("pop1", &["i1", "i2", "i1"])],
            &the_five_individuals(),
        )
        .unwrap_err();
        assert!(
            matches!(&error, Error::IndividualNamedTwiceInAPop { pop, name }
                if pop == "pop1" && name == "i1"),
            "{error:?}"
        );
    }

    /// A population with no individual has no value for any statistic,
    /// where pyNei gives NaN, so the name of that population is the error.
    #[test]
    fn a_population_with_no_individual_is_refused_with_its_name() {
        let error = Pops::from_names(
            &[pop_of("pop1", &["i1"]), pop_of("empty", &[])],
            &the_five_individuals(),
        )
        .unwrap_err();
        assert!(
            matches!(&error, Error::PopWithNoIndividual { pop } if pop == "empty"),
            "{error:?}"
        );
    }

    /// `pops` with no population at all would leave a result with nothing
    /// in it, where pyNei gives one with no column.
    #[test]
    fn pops_with_no_population_is_refused() {
        let error = Pops::from_names(&[], &the_five_individuals()).unwrap_err();
        assert!(matches!(&error, Error::NoPop), "{error:?}");
    }

    /// An individual that is in two populations is taken, as in pyNei,
    /// whose `test_maf_stats` names one individual in both of its
    /// populations. Each population holds it once.
    #[test]
    fn an_individual_in_two_populations_is_taken() {
        let pops = Pops::from_names(
            &[pop_of("pop1", &["i1", "i3"]), pop_of("pop2", &["i3", "i4"])],
            &the_five_individuals(),
        )
        .unwrap();
        assert_eq!(pops.individuals(0), [0, 2]);
        assert_eq!(pops.individuals(1), [2, 3]);
    }

    /// What a calculation that was given no `pops` works over: one
    /// population of every individual of the reader, in its order, named
    /// as pyNei names it.
    #[test]
    fn all_gives_one_population_named_pop_of_every_individual() {
        let pops = Pops::all(5);
        assert_eq!(pops.len(), 1);
        assert!(!pops.is_empty());
        assert_eq!(pops.name(0), DEFAULT_POP_NAME);
        assert_eq!(pops.name(0), "pop");
        assert_eq!(pops.individuals(0), [0, 1, 2, 3, 4]);
        assert!(pops.is_all(0));
    }

    /// A population a user named that turns out to be every individual in
    /// the order of the source is every individual: what `is_all` tells a
    /// caller is what the indices are, and not which constructor made
    /// them. One in another order, or one that leaves an individual out,
    /// is not.
    #[test]
    fn a_population_of_every_individual_in_the_order_of_the_source_is_all_of_them() {
        let pops = Pops::from_names(
            &[
                pop_of("in order", &["i1", "i2", "i3", "i4", "i5"]),
                pop_of("another order", &["i1", "i2", "i3", "i5", "i4"]),
                pop_of("four of them", &["i1", "i2", "i3", "i4"]),
            ],
            &the_five_individuals(),
        )
        .unwrap();
        assert!(pops.is_all(0));
        assert!(!pops.is_all(1));
        assert!(!pops.is_all(2));
    }
}
