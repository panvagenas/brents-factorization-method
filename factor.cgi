#!/usr/bin/perl
# Simple program that returns the prime factors, of a 
# given number, as json array.Input is readen from
# cgi param "number".
# Vagenas Panagiotis <pan.vagenas@gmail.com>
# https://bitbucket.org/vagenas/brents-factorization-method

# use strict;
# need CGI to fetch GET params
use CGI;
use Math::BigInt;

my $cgi = CGI->new();

printJsonOutput(factorize($cgi->param('number')));

exit;

# Finds factors of a given number
# input : any integer above 1
# returns : array containing the prime factors of input
sub factorize {
    my $input = Math::BigInt->new(@_);
    my @outPut;
    
    # check for valid input
    if ($input->is_nan() || $input->is_one() || $input->is_zero() || $input->is_neg()){
        return;
    }
    # check if input is a prime number and return it if so
    if(isPrime($input)){
        return $input;
    }
    # find a divisor of input
    my $divisor = Math::BigInt->new(brent($input));
    # recursively check divisor
    push(@outPut, factorize($divisor));
    # perform devision and recursively check result
    push(@outPut, factorize(Math::BigInt->new($input->bdiv($divisor))));
    
    @outPut;
}

# Richard's Brent factorization algorithm
# http://gan.anu.edu.au/~brent/pd/rpb051i.pdf
sub brent {
    my $N = Math::BigInt->new(@_);

    if ($N->copy()->bmod(2) == 0){
        return 2;;
    }
    # need POSIX to ceil random int
    use POSIX;
    my $c = Math::BigInt->new(ceil(rand($N->copy()->bsub(1))));
    my $m = Math::BigInt->new(ceil(rand($N->copy()->bsub(1))));
    my $y = Math::BigInt->new(ceil(rand($N->copy()->bsub(1))));

    my ($G, $r, $q) = (Math::BigInt->new, Math::BigInt->new(1), Math::BigInt->new(1));
    my ($k, $ys, $x, $min);

    do {
        $x=$y->copy();
        for (my $i = 0; $i < $r; $i++){
            $y->bmodpow(2,$N)->badd($c)->bmod($N);
        }
        $k = Math::BigInt->new(0);
        do {
            $ys = $y->copy();
            $min = $m->copy()->bsub($r)->is_pos() ? $r : $m;
            for (my $i = 0; $i < $min; $i++){
                $y->bmodpow(2,$N)->badd($c)->bmod($N);
                $q->bmul($x->copy()->bsub($y)->babs())->bmod($N);
            }
            $G = Math::BigInt::bgcd(($q, $N));
            $k->badd($m);
        } until ($k->bcmp($r) >= 0 || $G->bcmp(1) == 1);
        $r->bmul(2);
    } until ($G->bcmp(1) > 0);
    if (!$G->bcmp($N)){
        do {
            $ys->bmodpow(2,$N)->badd($c)->bmod($N);
            $G = Math::BigInt::bgcd(($x->copy()->bsub($ys)->babs(), $N));
        } until ($G->bcmp(1) > 0);
    }
    return $G;
}

# Checks if a given number is prime. First tries to use the AKS
# algorithm if Math::Primality::AKS module is installed, if not
# tries to load Math::Primality and makes use of the provided
# is_prime function that uses BPSW algorithim, at last if both
# attemts failed uses the implemented isMillerPrime function
# (we can't fully trust Miller–Rabin primality test until the
# generalized Riemann hypothesis is proved)
sub isPrime {
    use Module::Load::Conditional qw[can_load];
    my $n = Math::BigInt->new(@_);
    if (can_load( modules => {'Math::Primality::AKS' => undef})) {
        return is_aks_prime($n);
    } elsif (can_load( modules => {'Math::Primality' => undef})) {
        return is_prime($n);
    } else {
        return isMillerPrime($n, 60);
    }
}

# Miller-Rabin primality test
# http://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Algorithm_and_running_time
# www.cs.cornell.edu/courses/cs4820/2010sp/handouts/MillerRabin.pdf
sub isMillerPrime {
    my $n = Math::BigInt->new($_[0]);
    my $k = $_[1];

    return 1 if !$n->bcmp(2);
    return 0 if $n->bcmp(2) < 0 or !$n->copy()->bmod(2);
 
    my $d = $n->copy()->bsub(1);
    my $s = 0;
 
    while(!$d->copy()->bmod(2)) {
        $d->bdiv(2);
        $s++;
    }
    my $x;
    LOOP: for(1..$k) {
        $a = Math::BigInt->new(2 + int(rand($n->copy()->bsub(2))));
 
        $x = $a->bmodpow($d, $n);
        next if $x->is_one() or !$x->bcmp($n->copy()->bsub(1));
 
        for(1..$s-1) {
            $x->bmodpow(2,$n);
            return 0 if $x->is_one();
            next LOOP if !$x->bcmp($n->copy()->bsub(1));
        }
        return 0;
    }
    return 1;
}

# print json header and the given array as a json array
sub printJsonOutput {
    my (@out) = @_;
    my $response = join(',', @out);
    print $cgi->header(-type => "application/json", -charset => "utf-8");
    print qq{[$response]};
}