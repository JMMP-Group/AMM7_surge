#!/usr/bin/perl

use strict;
use warnings;

use Carp qw{croak};
use Fcm::Keyword::Config;
use Test::More (tests => 227);

BEGIN: {
    use_ok('Fcm::Keyword');
}

if (!caller()) {
    main(@ARGV);
}

sub main {
    local @ARGV = @_;
    local %Fcm::Keyword::Config::CONFIG_OF = (
        LOCATION_ENTRIES => {entry_class => 'Fcm::Keyword::Entry::Location'},
        REVISION_ENTRIES => {entry_class => 'Fcm::Keyword::Entry'},
    );
    test_get_prefix_of_location_keyword();
    test_get_entries();
    test_expand();
    test_unexpand();
    test_get_browser_url();
}

################################################################################
# Tests get_prefix_of_location_keyword().
sub test_get_prefix_of_location_keyword {
    is(Fcm::Keyword::get_prefix_of_location_keyword(), 'fcm');
    is(Fcm::Keyword::get_prefix_of_location_keyword(1), 'fcm:');
}

################################################################################
# Tests get_entries().
sub test_get_entries {
    my $entries = Fcm::Keyword::get_entries();
    isa_ok($entries, 'Fcm::Keyword::Entries');
    for (1 .. 10) {
        is(Fcm::Keyword::get_entries(), $entries, "get_entries: is singleton");
    }
    isnt(Fcm::Keyword::get_entries(1), $entries, "get_entries: can reset");
}

################################################################################
# Tests expand().
sub test_expand {
    my $T = 'expand';

    # Add some keywords for testing
    _add_keyword_entries([
    #   ['name'  , 'value'                   , {'rev1' => rev, ...}],
        ['FOO'   , 'test://foo/foo'      , {'V1.0' => 256, 'V1-1' => 4790}],
        ['FOO-TR', 'test://foo/foo/trunk', {}],
    ]);

    _do_keyword_tests($T, \&Fcm::Keyword::expand, [
    # Tests to ensure that valid targets are expanded
    #   [['input'               ], ['expected'                     ]],
        [['fcm:FOO'             ], ['test://foo/foo'               ]],
        [['fcm:FOO'     , 'V1.0'], ['test://foo/foo'       , '256' ]],
        [['fcm:Foo'             ], ['test://foo/foo'               ]],
        [['fcm:foo'             ], ['test://foo/foo'               ]],
        [['fcm:foo'     , 'v1.0'], ['test://foo/foo'       , '256' ]],
        [['fcm:foo'     , 'head'], ['test://foo/foo'       , 'head']],
        [['fcm:foo/'            ], ['test://foo/foo/'              ]],
        [['fcm:foo/'    , '1234'], ['test://foo/foo/'      , '1234']],
        [['fcm:foo/'    , 'v1.0'], ['test://foo/foo/'      , '256' ]],
        [['fcm:foo/'    , 'v1-1'], ['test://foo/foo/'      , '4790']],
        [['fcm:foo/bar'         ], ['test://foo/foo/bar'           ]],
        [['fcm:foo/bar' , 'PREV'], ['test://foo/foo/bar'   , 'PREV']],
        [['fcm:foo/bar' , 'base'], ['test://foo/foo/bar'   , 'base']],
        [['fcm:foo/bar' , 'v1-1'], ['test://foo/foo/bar'   , '4790']],
        [['fcm:foo/bar/', '7777'], ['test://foo/foo/bar/'  , '7777']],
        [['fcm:foo/bar/', '{11}'], ['test://foo/foo/bar/'  , '{11}']],
        [['fcm:foo/bar/', 'v1.0'], ['test://foo/foo/bar/'  , '256' ]],
        [['fcm:foo-tr'          ], ['test://foo/foo/trunk'         ]],
        [['fcm:foo-tr'  , 'head'], ['test://foo/foo/trunk' , 'head']],
        [['fcm:foo-tr'  , 'v1.0'], ['test://foo/foo/trunk' , '256' ]],
        [['fcm:foo-tr/'         ], ['test://foo/foo/trunk/'        ]],
        [['fcm:foo-tr/' , '1234'], ['test://foo/foo/trunk/', '1234']],
        [['fcm:foo-tr/' , 'v1-1'], ['test://foo/foo/trunk/', '4790']],
    # Tests to ensure that non-keyword targets are not expanded
    #   [['input'                     ]], # 'expected' same as 'input'
        [['no-change'                 ]],
        [['foo/bar'                   ]],
        [['/foo/bar'                  ]],
        [['/foo/bar'     , 'head'     ]],
        [['/foo/bar/'                 ]],
        [['/foo/bar/'    , 'not-a-key']],
        [['svn://foo/bar'             ]],
        [['svn://foo/bar', '1234'     ]],
        [['file://foo/bar'            ]],
        [['http://foo/bar'            ]],
    ]);

    # Tests for unexpected keywords
    for my $key (qw{foo bar baz}) {
        eval {
            Fcm::Keyword::expand("fcm:foo\@$key");
        };
        isa_ok($@, 'Fcm::Keyword::Exception', "$T: $key: invalid revision");
    }

    # Tests for "undef", all expecting exceptions
    for my $target_ref ([undef], [undef, undef], [undef, 'foo']) {
        eval {
            Fcm::Keyword::expand(@{$target_ref});
        };
        isa_ok($@, 'Fcm::Exception', "$T: undef");
    }
}

################################################################################
# Tests unexpand().
sub test_unexpand {
    my $T = 'unexpand';

    # Add some keywords for testing
    _add_keyword_entries([
    #   ['name'  , 'value'               , {'rev1' => rev, ...}],
        ['FOO'   , 'test://foo/foo'      , {'V1.0' => 256, 'V1-1' => 4790}],
        ['FOO_TR', 'test://foo/foo/trunk', {}],
        ['FOO-TR', 'test://foo/foo/trunk', {}],
    ]);

    _do_keyword_tests($T, \&Fcm::Keyword::unexpand, [
    # Tests to ensure that valid targets are expanded
    #   [['input'                        ], ['expected'            ]],
        [['test://foo/foo'               ], ['fcm:FOO'             ]],
        [['test://foo/foo'       , '256' ], ['fcm:FOO'     , 'V1.0']],
        [['test://foo/foo'       , 'head'], ['fcm:FOO'     , 'head']],
        [['test://foo/foo/'              ], ['fcm:FOO/'            ]],
        [['test://foo/foo/'      , '1234'], ['fcm:FOO/'    , '1234']],
        [['test://foo/foo/'      , '256' ], ['fcm:FOO/'    , 'V1.0']],
        [['test://foo/foo/'      , '4790'], ['fcm:FOO/'    , 'V1-1']],
        [['test://foo/foo/bar'           ], ['fcm:FOO/bar'         ]],
        [['test://foo/foo/bar'   , 'PREV'], ['fcm:FOO/bar' , 'PREV']],
        [['test://foo/foo/bar'   , 'base'], ['fcm:FOO/bar' , 'base']],
        [['test://foo/foo/bar'   , '4790'], ['fcm:FOO/bar' , 'V1-1']],
        [['test://foo/foo/bar/'  , '7777'], ['fcm:FOO/bar/', '7777']],
        [['test://foo/foo/bar/'  , '{11}'], ['fcm:FOO/bar/', '{11}']],
        [['test://foo/foo/bar/'  , '256' ], ['fcm:FOO/bar/', 'V1.0']],
        [['test://foo/foo/trunk'         ], ['fcm:FOO-TR'          ]],
        [['test://foo/foo/trunk' , 'head'], ['fcm:FOO-TR'  , 'head']],
        [['test://foo/foo/trunk' , '256' ], ['fcm:FOO-TR'  , 'V1.0']],
        [['test://foo/foo/trunk/'        ], ['fcm:FOO-TR/'         ]],
        [['test://foo/foo/trunk/', '1234'], ['fcm:FOO-TR/' , '1234']],
        [['test://foo/foo/trunk/', '4790'], ['fcm:FOO-TR/' , 'V1-1']],
    # Tests to ensure that non-keyword targets are not expanded
    #   [['input'                     ]], # 'expected' same as 'input'
        [['no-change'                 ]],
        [['foo/bar'                   ]],
        [['/foo/bar'                  ]],
        [['/foo/bar'     , 'head'     ]],
        [['/foo/bar/'                 ]],
        [['/foo/bar/'    , 'not-a-key']],
        [['svn://foo/bar'             ]],
        [['svn://foo/bar', '1234'     ]],
        [['file://foo/bar'            ]],
        [['http://foo/bar'            ]],
    ]);

    # Tests for "undef", all expecting exceptions
    for my $target_ref ([undef], [undef, undef], [undef, 'foo']) {
        eval {
            Fcm::Keyword::unexpand(@{$target_ref});
        };
        isa_ok($@, 'Fcm::Exception', "$T: undef");
    }
}

################################################################################
# Tests get_browser_url().
sub test_get_browser_url {
    my $T = 'get_browser_url';

    # Add some keywords for testing
    _add_keyword_entries([
    #   ['name'  , 'value'                        , {'rev1' => rev, ...}],
        ['FOO'   , 'test://foo/foo_svn/foo'       , {'V1' => 256, 'W2' => 479}],
        ['FOO-TR', 'test://foo/foo_svn/foo/trunk'],
        ['FOO_TR', 'test://foo/foo_svn/foo/trunk'],
    ]);

    my ($INPUT, $EXPECTED) = (0, 1);
    my ($LOC, $REV) = (0, 1);
    for my $test_ref (
    # Tests to ensure that valid targets are expanded
    #   [['input'                                ], 'expected'                                                ],
        [['test://foo/foo_svn/foo'               ], 'http://foo/projects/foo/intertrac/source:foo'            ],
        [['test://foo/foo_svn/foo'       , '256' ], 'http://foo/projects/foo/intertrac/source:foo@256'        ],
        [['test://foo/foo_svn/foo'       , 'head'], 'http://foo/projects/foo/intertrac/source:foo@head'       ],
        [['test://foo/foo_svn/foo/'              ], 'http://foo/projects/foo/intertrac/source:foo/'           ],
        [['test://foo/foo_svn/foo/'      , '1234'], 'http://foo/projects/foo/intertrac/source:foo/@1234'      ],
        [['test://foo/foo_svn/foo/'      , '256' ], 'http://foo/projects/foo/intertrac/source:foo/@256'       ],
        [['test://foo/foo_svn/foo/'      , '479' ], 'http://foo/projects/foo/intertrac/source:foo/@479'       ],
        [['test://foo/foo_svn/foo/bar'           ], 'http://foo/projects/foo/intertrac/source:foo/bar'        ],
        [['test://foo/foo_svn/foo/bar'   , '479' ], 'http://foo/projects/foo/intertrac/source:foo/bar@479'    ],
        [['test://foo/foo_svn/foo/bar/'  , '7777'], 'http://foo/projects/foo/intertrac/source:foo/bar/@7777'  ],
        [['test://foo/foo_svn/foo/bar/'  , '{11}'], 'http://foo/projects/foo/intertrac/source:foo/bar/@{11}'  ],
        [['test://foo/foo_svn/foo/bar/'  , '256' ], 'http://foo/projects/foo/intertrac/source:foo/bar/@256'   ],
        [['test://foo/foo_svn/foo/trunk'         ], 'http://foo/projects/foo/intertrac/source:foo/trunk'      ],
        [['test://foo/foo_svn/foo/trunk' , 'head'], 'http://foo/projects/foo/intertrac/source:foo/trunk@head' ],
        [['test://foo/foo_svn/foo/trunk' , '256' ], 'http://foo/projects/foo/intertrac/source:foo/trunk@256'  ],
        [['test://foo/foo_svn/foo/trunk/'        ], 'http://foo/projects/foo/intertrac/source:foo/trunk/'     ],
        [['test://foo/foo_svn/foo/trunk/', '1234'], 'http://foo/projects/foo/intertrac/source:foo/trunk/@1234'],
        [['test://foo/foo_svn/foo/trunk/', '479' ], 'http://foo/projects/foo/intertrac/source:foo/trunk/@479' ],
        [['fcm:FOO'                              ], 'http://foo/projects/foo/intertrac/source:foo'            ],
        [['fcm:FOO'                      , 'V1'  ], 'http://foo/projects/foo/intertrac/source:foo@256'        ],
        [['fcm:FOO'                      , 'head'], 'http://foo/projects/foo/intertrac/source:foo@head'       ],
        [['fcm:FOO/'                             ], 'http://foo/projects/foo/intertrac/source:foo/'           ],
        [['fcm:FOO/'                     , '1234'], 'http://foo/projects/foo/intertrac/source:foo/@1234'      ],
        [['fcm:FOO/'                     , 'V1'  ], 'http://foo/projects/foo/intertrac/source:foo/@256'       ],
        [['fcm:FOO/'                     , 'W2'  ], 'http://foo/projects/foo/intertrac/source:foo/@479'       ],
        [['fcm:FOO/bar'                          ], 'http://foo/projects/foo/intertrac/source:foo/bar'        ],
        [['fcm:FOO/bar'                  , 'W2'  ], 'http://foo/projects/foo/intertrac/source:foo/bar@479'    ],
        [['fcm:FOO/bar/'                 , '7777'], 'http://foo/projects/foo/intertrac/source:foo/bar/@7777'  ],
        [['fcm:FOO/bar/'                 , '{11}'], 'http://foo/projects/foo/intertrac/source:foo/bar/@{11}'  ],
        [['fcm:FOO/bar/'                 , 'v1'  ], 'http://foo/projects/foo/intertrac/source:foo/bar/@256'   ],
        [['fcm:FOO-TR'                           ], 'http://foo/projects/foo/intertrac/source:foo/trunk'      ],
        [['fcm:FOO-TR'                   , 'head'], 'http://foo/projects/foo/intertrac/source:foo/trunk@head' ],
        [['fcm:FOO-TR'                   , 'V1'  ], 'http://foo/projects/foo/intertrac/source:foo/trunk@256'  ],
        [['fcm:FOO-TR/'                          ], 'http://foo/projects/foo/intertrac/source:foo/trunk/'     ],
        [['fcm:FOO-TR/'                  , '1234'], 'http://foo/projects/foo/intertrac/source:foo/trunk/@1234'],
        [['fcm:FOO-TR/'                  , 'w2'  ], 'http://foo/projects/foo/intertrac/source:foo/trunk/@479' ],
    ) {
        my $input = $test_ref->[$INPUT][$LOC];
        if (exists($test_ref->[$INPUT][$REV])) {
            $input .= '@' . $test_ref->[$INPUT][$REV];
        }
        for (
            {name => "$T: scalar input: $input", input => [$input]},
            {name => "$T: list input: $input"  , input => $test_ref->[$INPUT]},
        ) {
            my $output;
            eval {
                $output = Fcm::Keyword::get_browser_url(@{$_->{input}});
                is($output, $test_ref->[$EXPECTED], $_->{name});
            };
            if ($@) {
                fail("$_->{name}: $@");
            }
        }
    }

    # Tests correct behaviour for "undef"
    for my $bad_url (undef, '') {
        eval {
            Fcm::Keyword::get_browser_url($bad_url);
        };
        isa_ok($@, 'Fcm::Exception', sprintf(
            "$T: %s", (defined($bad_url) ? $bad_url : 'undef'),
        ));
    }

    # Tests correct behaviour for invalid inputs
    for my $bad_url ('foo', 'svn://no/such/url', 'fcm:no_such_project/trunk') {
        eval {
            Fcm::Keyword::get_browser_url($bad_url);
        };
        isa_ok($@, 'Fcm::Keyword::Exception', "$T: $bad_url: invalid keyword");
    }
}

################################################################################
# Adds keyword entries.
sub _add_keyword_entries {
    my ($items_ref) = @_;
    my ($NAME, $LOC, $REV) = (0 .. 2);
    my $entries = Fcm::Keyword::get_entries(1); # reset
    for my $item_ref (@{$items_ref}) {
        my $entry = $entries->add_entry($item_ref->[$NAME], $item_ref->[$LOC]);
        while (my ($key, $value) = each(%{$item_ref->[$REV]})) {
            $entry->get_revision_entries()->add_entry($key, $value);
        }
    }
}

################################################################################
# Performs keyword testings.
sub _do_keyword_tests {
    my ($T, $action_ref, $tests_ref) = @_;
    my ($INPUT, $EXPECTED) = (0, 1);
    my ($LOC, $REV) = (0, 1);
    for my $test_ref (@{$tests_ref}) {
        if (!defined($test_ref->[$EXPECTED])) {
            $test_ref->[$EXPECTED] = $test_ref->[$INPUT];
        }
        my %value_of;
        for my $i (0 .. $#{$test_ref}) {
            $value_of{$i} = $test_ref->[$i][$LOC];
            if (exists($test_ref->[$i][$REV])) {
                $value_of{$i} .= '@' . $test_ref->[$i][$REV];
            }
        }
        eval {
            is(
                $action_ref->($value_of{$INPUT}), $value_of{$EXPECTED},
                "$T: scalar context: $value_of{$INPUT}",
            );
        };
        if ($@) {
            fail("$T: scalar context: $value_of{$INPUT}: $@");
        }
        eval {
            is_deeply(
                [$action_ref->(@{$test_ref->[$INPUT]})],
                $test_ref->[$EXPECTED],
                "$T: list context: $value_of{$INPUT}",
            );
        };
        if ($@) {
            fail("$T: list context: $value_of{$INPUT}: $@");
        }
    }
}

__END__
