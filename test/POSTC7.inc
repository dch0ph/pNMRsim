POSTC7elem {
	pulse $1 $2 $3
	pulse ($1)*4 $2 ($3)+180
	pulse ($1)*3 $2 $3
}

C7elem {
	pulse ($1)*8 $2 $3
}

C7 {
    include C7elem $1 $2 $3+0
    include C7elem $1 $2 ($3)+51.4286
   include C7elem $1 $2 ($3)+102.8571
   include C7elem $1 $2 ($3)+154.2857
   include C7elem $1 $2 ($3)+205.7143
    include C7elem $1 $2 ($3)+257.1429
    include C7elem $1 $2 ($3)+308.5714
}

POSTC7 {
    include POSTC7elem $1 $2 $3+0
    include POSTC7elem $1 $2 ($3)+51.4286
   include POSTC7elem $1 $2 ($3)+102.8571
   include POSTC7elem $1 $2 ($3)+154.2857
   include POSTC7elem $1 $2 ($3)+205.7143
    include POSTC7elem $1 $2 ($3)+257.1429
    include POSTC7elem $1 $2 ($3)+308.5714
}
