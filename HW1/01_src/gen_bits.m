function b = gen_bits(seed)
% Generate the 32 random bits b_0 ~ b_31 shared by all HW1 questions
rng(seed, 'twister');
b = randi([0 1], 1, 32);
end
