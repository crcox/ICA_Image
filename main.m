image_list = {
    'airplane_bw_small.png'
    'arrow_bw_small.png'
    'bell_bw_small.png'
    'book_bw_small.png'
    'bottle_bw_small.png'
    'cat_bw_small.png'
    'chicken_bw_small.png'
    'cow_bw_small.png'
    'crow_bw_small.png'
    'cuckoo_clock_bw_small.png'
    'dog_bw_small.png'
    'drum_bw_small.png'
    'elephant_bw_small.png'
    'explosion_bw_small.png'
    'fire_bw_small.png'
    'frog_bw_small.png'
    'glass_bw_small.png'
    'goat_bw_small.png'
    'goose_bw_small.png'
    'guitar_bw_small.png'
    'gun_bw_small.png'
    'helicopter_bw_small.png'
    'horse_bw_small.png'
    'key_bw_small.png'
    'lion_bw_small.png'
    'mobilephone_bw_small.png'
    'motorcycle_bw_small.png'
    'piano_bw_small.png'
    'rain_bw_small.png'
    'scissors_bw_small.png'
    'sword_bw_small.png'
    'tennis_bw_small.png'
    'toilet_bw_small.png'
    'train_bw_small.png'
    'whistle_bw_small.png'
    'wolf_bw_small.png'
    'zipper_bw_small.png'
};

%%
imagedir = 'D:\ECoG\KyotoNaming\Rick\stimuli\Pictures_Snodgrass';
x = dir(imagedir);
animate = [0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,1,0,1,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1];
[~,anisort] = sort(animate);
z = ~cellfun('isempty', regexp({x.name},'.*.bmp$'));
image_list = {x(z).name};

IMAGES = cell(numel(image_list),1);
for i = 1:numel(IMAGES)
    IMAGES{i} = imresize(imread(fullfile(imagedir,image_list{i}), 'bmp'),[150,150]);
end


%%
nw = 5;
dim = 150;
windows = cell(nw);
[a,b] = ndgrid(0:(nw-1),0:(nw-1));
a = (a(:) .* (150/nw));
b = (b(:) .* (150/nw));
n = 150 / nw;
for i = 1:numel(windows)
    w = false(150);
    [aix,bix] = ndgrid(a(i) + (1:n), b(i) + (1:n));
    w(aix(:),bix(:)) = true;
    windows{i} = w;
end

getWindow = @(y) cell2mat(cellfun(@(x) x(windows{y}), IMAGES, 'UniformOutput', false)');

%%
% addpath('/Users/Chris/Documents/MATLAB/Toolboxes/FastICA_25');
addpath('C:\Users\mbmhscc4\MATLAB\Toolboxes\FastICA_25');
X = zeros(900, 100 * numel(windows));
cur = 0;
for i = 1:numel(windows)
    disp(i)
    a = cur + 1;
    b = cur + 100;
    cur = b;
    X(:,a:b) = getWindow(i);
%     disp(tabulate(X(:)))
end
[IC,A,W] = fastica(X,'numOfIC',900);
IC = mat2cell(IC, 900, repmat(100,1,numel(windows)));

%%
getImg = @(X,y) cell2mat(cellfun(@(x) reshape(x(:,y),n,n), X, 'UniformOutput', 0));
IMAGE_RECON = cell(numel(image_list),1);
X = cell(numel(windows),1);
for i = 1:numel(windows)
    w = windows{i};
    X{i} = A*IC{i};
end
% X{10} = zeros(size(X{1}));
for j = 1:numel(IMAGE_RECON)
    IMAGE_RECON{j} = getImg(reshape(X,nw,nw),j);
end
for i = 1:36, subplot(6,6,i); imagesc(IMAGE_RECON{i}); end

%%
U = cell(nw);
S = cell(nw);
V = cell(nw);
for i = 1:numel(U)
    disp(i);
    X = double(getWindow(i));
%     disp(tabulate(X(:)))
    try
        [U{i},S{i},V{i}] = svd(X);
    catch
        % nothing
    end
end

%%
getImg = @(X,y) cell2mat(cellfun(@(x) reshape(x(:,y),n,n), X, 'UniformOutput', 0));
IMAGE_RECON = cell(numel(image_list),1);
X = cell(numel(windows),1);
d = 10;
for i = 1:numel(windows)
    w = windows{i};
    s = diag(S{i}(1:numel(IMAGES),:));
    s((d+1):end) = 0;
    Sd = S{i};
    Sd(1:numel(IMAGES),:) = diag(s);
    X{i} = U{i}*Sd*V{i}';
end
for j = 1:numel(IMAGE_RECON)
    IMAGE_RECON{j} = getImg(reshape(X,nw,nw),j);
end
for i = 1:36, subplot(6,6,i); imagesc(IMAGE_RECON{i}); end

imagesc(squareform(pdist(tmp(anisort,:))))